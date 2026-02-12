### Friday December 26th, 2025
# This is an R script which takes as input file(s) to shuffle over TE annotations to calculate enrichments for each TE subfamily.
# The output will be used to make all the cool enrichment plots using mm10. 
# 
# Expected script call:
# Rscript --vanilla TEshuffle.r [file input directory] [file output directory] [pattern to match] [cores]

# Load in all the packages you need
suppressPackageStartupMessages({ 
	library(data.table); library(ChIPseeker); library(doParallel); library(regioneR); library(tidyr); library(dplyr); 
	library(GenomicFeatures)
})

# These are the cool functions you made to generate the ability to overlap and shuffle the intervals, respectively. Yeah!
overlap <- function(A, B, minFrac, aBool, bBool, ...) {
	if(!hasArg(minFrac)) minFrac <- 1; if(!hasArg(aBool)) aBool <- T; if(!hasArg(bBool)) bBool <- T

	# Hits can return multiple intersections, so if an entry in A is intersected twice by B, it will have two rows - one for each intersection
	hits <- findOverlaps(A, B)

	# Overlaps returns a subset of A which has at least one intersection with B. 
	overlaps <- pintersect(A[queryHits(hits)], B[subjectHits(hits)])

	if (aBool & bBool) { 
		return( queryHits(hits[(width(overlaps) / width(B[subjectHits(hits)])) >= minFrac | (width(overlaps) / width(A[queryHits(hits)])) >= minFrac]) )
	} else if (aBool) { 
		return( queryHits(hits[(width(overlaps) / width(A[queryHits(hits)])) >= minFrac]) ) 
	} else { return( queryHits(hits[(width(overlaps) / width(B[subjectHits(hits)])) >= minFrac]) ) }
}

shuffle <- function(A, chr_sizes, tss, tss_chr) {
	
	output <- GRanges()
	a_chr <- as.data.frame(table(seqnames(A))); rownames(a_chr) <- a_chr$Var1
	for (f in 1:nrow(a_chr)) {
		# Get the current chromosome name to iterate through
		chr.name <- as.character(a_chr[f,1])

		# Subset out those elements in the current family that are on the same chromosome
		test <- A[ which(seqnames(A) == chr.name), ]

		# Now generate which random TSS each will be shuffled to and create a subset
		rand.tss <- round(runif(length(test), 1, tss_chr[ chr.name, 2 ])); tss.sub <- tss[ which(seqnames(tss) == chr.name), ][ rand.tss, ]

		# Generate the new shuffled coordinates for each interval based on the random TSS and their orrientation
		# using the relative strand: + (left of TSS) or - (right of TSS) 
		# NOTE: Since this is summit coords, end - start yields 0, but when *not* using summits this value will be nonzero 
		test$newStart <- ifelse(test$rel_strand == "+", start(tss.sub) - test$distance - (end(test) - start(test)), end(tss.sub) + test$distance) 
		test$newEnd <- test$newStart + (end(test) - start(test))

		# Check whether any intervals are above or below the chromosome coordinates 
		vals <- which(test$newStart <= 1 | test$newEnd >= (chr_sizes[ chr.name, 1 ]-1))
	
		# If there are intervals that need changed, change them according to their relative strand 
		if (length(vals) != 0) { 
			tmp <- test[vals, ]
			test[vals,]$newStart <- ifelse(tmp$rel_strand == "+", end(tss.sub)[vals] + tmp$distance, start(tss.sub)[vals] - tmp$distance - (end(tmp) - start(tmp)))  
			test[vals,]$newEnd <- test[vals,]$newStart + (end(test[vals,]) - start(test[vals,])) }

		output <- c(output, test)
	}
	ranges(output) <- IRanges(start = output$newStart, end = output$newEnd); return( output )
}

# Load in the arguements, such as the current path you want to load the files from, the path you want to output results to, and the number of
# cores you want to run
args <- commandArgs(TRUE) 

# Assume all 4 arguements are provided
path <- args[1]; res_path <- args[2]; pat <- args[3]; cores.main <- args[4]
dir.create(res_path)

# Just keep this at 1 for now; there's no main reason to have multiple samples going simultaneously
cores.second <- 1

# Load in the current genome sizes and TSS sites. 
chrMax <- 19
sizes <- read.table("../STAR_noGene_mm10/chrNameLength.txt", row.names = 1); colnames(sizes)[1] <- "size"
txdb <- GenomicFeatures::makeTxDbFromGFF( list.files(pattern = "mm10.refseqCurated.gtf.*")[1] )
tss <- promoters(txdb, upstream = 0, downstream = 3)

# Generate a table of the total number of TSS sites on each chr
tss_chr <- as.data.frame(table(seqnames(tss))); rownames(tss_chr) <- tss_chr$Var1

# Load in the TE file
TEs <- read.table("rm/mm10_repeatMasker.txt.gz"); colnames(TEs)[1:3] <- c("chr", "start", "end")

# Get the data.frame ready to iterate through with the permTest
TE_families <- (unique(TEs$V4))

# Turn the TEs into a GRanges for faster analysis
TEs <- makeGRangesFromDataFrame(TEs, starts.in.df.are.0based=T, keep.extra.columns=T)
seqlevels(TEs) <- sort(seqlevels(TEs))

# Make a data.frame that can be used as a base for each analysis below
fam_base <- data.frame(group = TE_families, value = rep(0, length(TE_families)))

# Create your own famStats file, since that's easier than having to input one every time...
famStats <- as.data.frame(table(paste0(TEs$V4, "|", TEs$V5))); famStats <- separate(famStats, col = "Var1", into = c("repName", "repClass_Family"), sep = "\\|")
TEs <- subset(TEs, !(V4 %in% subset(famStats, Freq <= 20 | grepl("Satellite", famStats$repClass_Family) | grepl("Unknown", famStats$repClass_Family) )$repName))

# Use the path from above to collect a list of files that you'll be working on to loop over. 
files <- list.files(path = path, pattern = pat) 

TEs$id <- paste0(TEs$V4, "/", TEs$V5)
results <- mclapply(1:length(files), function(i) {

	# Make a name variable for saving to the file
	name <- sub( sub("^.", "", pat), "", files[i])

	# Load in the peaks you need and make sure they're all filtered and sorted
	peaks <- read.table(paste0(path, "/", files[i]), sep = "\t"); colnames(peaks)[grep(c(paste0("V", 1:3, "$", collapse = "|")), colnames(peaks))] <- c("chr", "start", "end")
	peaks <- makeGRangesFromDataFrame(peaks, starts.in.df.are.0based=T, keep.extra.columns=T); peaks <- filterChromosomes(peaks, keep.chr=unlist(c(paste0("chr", c(1:chrMax)), "chrX")))
	seqlevels(peaks) <- sort(seqlevels(peaks)); peaks <- sort(peaks)

	# Save the above output to a data.frame and then filter the annotations a bit to make for ease of writing to table.
	a <- as.data.frame(annotatePeak(peaks, tssRegion=c(-1000,1000), TxDb=txdb, verbose = F))
	a$anno <- a$annotation; a[grepl("Exon", a$anno), "anno"] <- "Exon"; a[grepl("Intron", a$anno), "anno"] <- "Intron"

	# Use the summit information called in MACS to only have the summits used for a simpler shuffle
	peaks <- GRanges(seqnames = as.character(seqnames(peaks)), ranges = IRanges(c(start(peaks) + peaks$V10), width = 1), strand = strand(peaks)) 

	# Append the annotation (in case summit is T and would erase the info)
	peaks$anno <- ifelse( a$anno == "Promoter", "Promoter", ifelse( a$anno %in% c("Distal Intergenic", "Downstream (<=300bp)", "Intron"), "Enhancer", "Gene-Sequence" ) )

	# Generate a vector of distances to the nearest TSS for the current peak file. Include the relative strand (rel_strand) info to try 
	# to put the particular peaks in the same relative position against the TSS when you shuffle them
	dist <- as.data.frame(distanceToNearest(peaks, tss)); dist$rel_strand <- ifelse(start(tss[dist$subjectHits,]) - end(peaks[dist$queryHits,]) >= 0, "+", "-")
	peaks$distance <- dist$distance; peaks$rel_strand <- dist$rel_strand
	
	# If there are <= 10 overlaps total for either enhancers or promoters, skip the current factor and move on to the next one
	df <- data.frame()
	sub <- subset(peaks, anno == "Enhancer")
	olap <- overlap(TEs, sub, 1)
	if (length(olap) <= 10) { 
		init_fams <- data.frame(repName = NA, repClass = NA, repFamily = NA, total = length(sub), bound_uniq = NA, bound = length(olap), exp_uniq = NA, exp = NA, perm_enriched_pval = NA, fisher_enriched_padj = NA, perm_depleted_pval = NA, fisher_depleted_padj = NA, samp = paste0(name, "_only", "Enhancer"))
		write.table(init_fams, file = paste0(res_path, "/", name, "_onlyEnhancer", ".txt"), row.names = F, quote = F, sep = "\t")
		df <- rbind(df, init_fams)
	} else {
		# Make a temp file which contains the initial family values for intersections against sub
		te.dat <- as.data.frame(TEs); te.dat$te.dat <- 0; te.dat$id <- paste0(te.dat$V4, "/", te.dat$V5)
		te.dat[unique(as.numeric(as.character(olap))),"te.dat"] <- as.data.frame(table(olap))[,2]

		# Now combine the above te.dat file into the init_fams file for family stat calculations later!
		init_fams <- group_by(te.dat, id) %>% summarise(total = n(), bound_uniq = sum(te.dat >= 1), bound = sum(te.dat)) %>% as.data.frame()
		init_fams$exp <- init_fams$exp_uniq <- 0;  init_fams$fisher_depleted_padj <- init_fams$perm_depleted_pval <- init_fams$fisher_enriched_padj <- init_fams$perm_enriched_pval <- 1

		# Subset out TEs based on the family size and number of bound elements (total number is <= 100 & no "bound" elements) 
		init_fams <- subset(init_fams, !(total <= 100 & bound == 0) | total > 100); subTEs <- subset(TEs, id %in% init_fams$id)
		te.dat <- subset(te.dat, id %in% init_fams$id)

		# The main shuffle method right here! Run for 1000 permutations with "cores.main" cores for each time
		output <- mclapply(1:1000, function(f) {
			res <- overlap(subTEs, shuffle(sub, sizes, tss, tss_chr), 1); tmp <- te.dat; tmp$tmp <- 0
			if (length(res) != 0) { tmp[unique(as.numeric(as.character(res))),"tmp"] <- as.data.frame(table(res))[,2] }
			group_by(tmp, id) %>% summarise(total = n(), bound_uniq = sum(tmp >= 1), bound = sum(tmp)) %>% as.data.frame()	
		}, mc.cores = cores.main)

		# Iterate over all of the different results from the above permuations and measure whether there is 2X more bound than expected
		for (f in 1:length(output)) {
			init_fams$exp_uniq <- init_fams$exp_uniq + output[[f]]$bound_uniq; init_fams$exp <- init_fams$exp + output[[f]]$bound
			init_fams$perm_enriched_pval <- ifelse(init_fams$bound_uniq > (2 * output[[f]]$bound_uniq), init_fams$perm_enriched_pval, init_fams$perm_enriched_pval + 1)
			init_fams$perm_depleted_pval <- ifelse(init_fams$bound_uniq < (output[[f]]$bound_uniq / 2), init_fams$perm_depleted_pval, init_fams$perm_depleted_pval + 1) 
		}

		# Now calculate the permuation pvals by dividing both by output length + 1
		init_fams$perm_enriched_pval <- init_fams$perm_enriched_pval / (length(output)+1); init_fams$perm_depleted_pval <- init_fams$perm_depleted_pval / (length(output)+1)

		# Divide expected by the number of permutations to get the mean expected number of sub within each interval
		init_fams$exp_uniq <- init_fams$exp_uniq / length(output); init_fams$exp <- init_fams$exp / length(output)

		# Calculate the fisher significance for enrichment and depletion and adjust the p-values using p.adjust
		init_fams$fisher_enriched_padj <- p.adjust(unlist(apply(init_fams, 1, function(x) fisher.test(matrix(c(as.numeric(x["bound_uniq"]), as.numeric(x["total"]) - as.numeric(x["bound_uniq"]), round(as.numeric(x["exp_uniq"])), round(as.numeric(x["total"]) - as.numeric(x["exp_uniq"]))), nrow = 2), alternative = "greater")$p.value)), method = "BH")
		init_fams$fisher_depleted_padj <- p.adjust(unlist(apply(init_fams, 1, function(x) fisher.test(matrix(c(as.numeric(x["bound_uniq"]), as.numeric(x["total"]) - as.numeric(x["bound_uniq"]), round(as.numeric(x["exp_uniq"])), round(as.numeric(x["total"]) - as.numeric(x["exp_uniq"]))), nrow = 2), alternative = "less")$p.value)), method = "BH")

		# Add in the sample name for later merging!
		init_fams$samp <- paste0(name, "_onlyEnhancer")

		# Separate out the "id" column for repName, repClass, and repFamily
		init_fams <- separate(init_fams, col = "id", sep = "\\/", into = c("repName", "repClass", "repFamily"))

		# Instead of returning the object, save it as an .txt object 
		write.table(arrange(init_fams, fisher_enriched_padj, desc(bound)), file = paste0(res_path, "/", name, "_onlyEnhancer", ".txt"), row.names = F, quote = F, sep = "\t") 

		# Now return the final object for each sample. Yay!
		df <- rbind(df, init_fams)
	}

	return(df)
}, mc.cores = cores.second)

# With the above complete, generate a results object and write it to disk so that you can access it later!
write.table(arrange(bind_rows(results), fisher_enriched_padj, desc(bound)), file = paste0(res_path, "/TEshuffle_results.txt"), sep = "\t", quote = F, row.names = F)
system( paste("gzip", paste0(res_path, "/TEshuffle_results.txt")) ) 




