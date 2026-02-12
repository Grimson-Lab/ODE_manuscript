### Friday December 26th, 2025
# This is a copy of your R TEshuffle script. Instead, it will be used for testing whether there's significant overlap between the 
# ABC-linkages and the ATAC-RNA correlation-derived linkages. Woot woot
# 
# Expected script call (post running module load R/4.1.3-r9):
# Rscript --vanilla TEshuffle.r [input file] [cores]

# Load in all the packages you need
suppressPackageStartupMessages({ 
	library(data.table); library(ChIPseeker); library(doParallel); library(regioneR); library(tidyr); library(dplyr)
 
})

# These are the cool functions you made to generate the ability to overlap and shuffle the intervals, respectively. Yeah!
overlap <- function(A, B, ...) {

	# Hits can return multiple intersections, so if an entry in A is intersected twice by B, it will have two rows - one for each intersection
	hits <- findOverlaps(A, B)

	# From the above, figure out which of these hits are "valid", i.e. the link is to the same TSS for both A and B!
	# This requires that the used TSS from "A" 1) exists in the "lut" object, and also 2) the "peak_id" from B is in the list of identified
	# links from "lut"
	tmp <- as.data.frame(hits); tmp$A_tss <- A[queryHits(hits), ]$tss
	tmp$inLut <- sapply(1:nrow(tmp), function(i) { nrow(lut[.(tmp[i,"A_tss"]), nomatch = 0]) > 0 } )
	tmp$B_peakID <- B[subjectHits(hits), ]$peak_id
	tmp$peakInLut <- sapply(1:nrow(tmp), function(i) { if (tmp[i,"inLut"]) { if ( tmp[i,"B_peakID"] %in% lut[J(tmp[i,"A_tss"]),val][[1]] ) { return( TRUE ) } else { return (FALSE) } } else { return(NA) } } )

	return( nrow( unique(subset(tmp, peakInLut == T)[,c("A_tss", "B_peakID")]) ) )
}

shuffle <- function(A, chr, chr_sizes, tss, tss_chr) {

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
		test$tss <- tss[ rand.tss, ]$id
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

# Assume all arguments are provided
input <- args[1]; cores.main <- args[2]

# Since the input needs to be edited a bit, do that below
colnames(input) <- c("chr", "start", "end", "peak_id", "tss_pos", "tss")
input$distance <- sapply(1:nrow(input), function(i) { min( abs( input[i,"tss_pos"] - input[i,"end"] ), abs( input[i,"tss_pos"] - input[i,"start"] ) ) })
input$rel_strand <- ifelse( input$tss_pos < input$start, "-", "+" )

input <- makeGRangesFromDataFrame(input, starts.in.df.are.0based=T, keep.extra.columns=T)
seqlevels(input) <- sort(seqlevels(input)); input <- sort(input)

# Load in the current genome sizes and TSS sites. 
sizes <- read.table("../STAR_noGene_mm10/chrNameLength.txt", row.names = 1); colnames(sizes)[1] <- "size"
txdb <- GenomicFeatures::makeTxDbFromGFF("../mm10.ncbiRefSeq.gtf.gz")
tss <- readRDS(file = "atac_rna_cor_tss.rds")

# Generate a table of the total number of TSS sites on each chr
tss_chr <- as.data.frame(table(seqnames(tss))); rownames(tss_chr) <- tss_chr$Var1

# Load in the look-up table object to check which peaks are linked with each TSS from the ATAC-RNA-correlation analysis
lut <- readRDS(file = "atac_rna_cor_lut.rds")

# Load in the intervals for the ATAC-RNA linkages and make sure they're sorted
arl <- read.table("atac_rna_linkages.txt")[,c(5:8,4)]; colnames(arl) <- c("chr", "start", "end", "peak_id", "gene")
arl <- makeGRangesFromDataFrame(arl, starts.in.df.are.0based=T, keep.extra.columns=T)
seqlevels(arl) <- sort(seqlevels(arl)); arl <- sort(arl)

# Determine the initial number of overlaps between the two datasets
olap <- overlap(input, arl)

# Now combine the above te.dat file into the init_fams file for family stat calculations later!
init <- data.frame("total" = length(input), "overlap" = olap, "exp" = 0)
init$fisher_depleted_padj <- init$perm_depleted_pval <- init$fisher_enriched_padj <- init$perm_enriched_pval <- 1

# The main shuffle method right here! Run for 1000 permutations with "cores.main" cores for each time
output <- mclapply(1:1000, function(f) {
	overlap( shuffle(input, T, sizes, tss, tss_chr), arl )	
}, mc.cores = cores.main)

# Iterate over all of the different results from the above permuations and measure whether there is 2X more overlap than expected
for (f in 1:length(output)) {
	init$exp <- init$exp + output[[f]]
	init$perm_enriched_pval <- ifelse(init$overlap > (2 * output[[f]]), init$perm_enriched_pval, init$perm_enriched_pval + 1)
	init$perm_depleted_pval <- ifelse(init$overlap < (output[[f]] / 2), init$perm_depleted_pval, init$perm_depleted_pval + 1) 
}

# Now calculate the permuation pvals by dividing both by output length + 1
init$perm_enriched_pval <- init$perm_enriched_pval / (length(output)+1); init$perm_depleted_pval <- init$perm_depleted_pval / (length(output)+1)

# Divide expected by the number of permutations to get the mean expected number of sub within each interval
init$exp <- init$exp / length(output)

# Calculate the fisher significance for enrichment and depletion and adjust the p-values using p.adjust
init$fisher_enriched_padj <- unlist(apply(init, 1, function(x) fisher.test(matrix(c(as.numeric(x["overlap"]), as.numeric(x["total"]) - as.numeric(x["overlap"]), round(as.numeric(x["exp"])), round(as.numeric(x["total"]) - as.numeric(x["exp"]))), nrow = 2), alternative = "greater")$p.value))
init$fisher_depleted_padj <- unlist(apply(init, 1, function(x) fisher.test(matrix(c(as.numeric(x["overlap"]), as.numeric(x["total"]) - as.numeric(x["overlap"]), round(as.numeric(x["exp"])), round(as.numeric(x["total"]) - as.numeric(x["exp"]))), nrow = 2), alternative = "less")$p.value))

# Instead of returning the object, save it as an .txt object 
write.table(init, file = "abc_enrichment_in_correlation_links.txt", row.names = F, quote = F, sep = "\t") 
