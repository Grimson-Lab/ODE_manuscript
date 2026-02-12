# This is a compendium of the R code used to generate the ~460 million correlations between ATAC-seq peaks and genes on the same chromosome
# for mm10 using the Immgen publicly available immune cell type data. It will take a *while* to do that, so save plenty of cores and time

### NOTE: Make sure to add "lbzip2" to your environment prior to starting R!

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(data.table); setDTthreads(8); library(dplyr); library(fastSave); library(GenomicRanges); library(doParallel)
})

# Begin by reading in the list of RefSeq Curated gene annotations for mm10. Yeah
d <- subset(read.table("ncbiRefSeqCurated.txt.gz", sep = "\t"), !grepl("chrM|_alt|chrUn|_random|_fix|chrY", V3) )[,c("V3", "V5", "V6", "V13", "V4")]
d$end <- d$start <- ifelse( d$V4 == "+", d$V5, d$V6 ); d <- unique(d[,c("V3", "start", "end", "V13")])
d <- as.data.table(d); setkey(d, V13)

# Turn the list of *ALL* peaks to permutate into a GRanges for faster analysis
open <- as.data.frame(fread("atac/summit_unif_peaks.bed")); colnames(open) <- c("chr", "start", "end", "id")
open <- makeGRangesFromDataFrame(open, starts.in.df.are.0based=T, keep.extra.columns=T); open$coords <- with(as.data.frame(open), paste0(seqnames, ":", start, "-", end))
seqlevels(open) <- sort(seqlevels(open)); open <- sort(open)

# Load in the quantile normalized ATAC-seq values
load.lbzip2(file = "atac/rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)

# Load in the RNA-seq output, and also do quantile normalization on it to get fancy comparable numbers. Woot woot?
gex <- read.csv("rna/refseqCurated_tpm.csv", header = T, row.names = 1)

# Generate informative objects for more rapidly performing correlation tests
l2.gex <- log2(gex[ , colnames(l2.atac.cpm.qnorm)[colnames(l2.atac.cpm.qnorm) %in% colnames(gex)]]+1)

# Temporarily write the above object so that it can be qnormed 
write.table(l2.gex, file = "rna/rObject_l2.gex.txt", sep = "\t", quote = F)

# Now, in unix, call "qnorm" to generate the quantile normalized values on the above object
system("qnorm rna/rObject_l2.gex.txt > rna/rObject_l2.gex.qnorm.txt")

# Load the qnorm output 
l2.gex.qnorm <- read.table("rna/rObject_l2.gex.qnorm.txt", header = T, row.names = 1)

# Save the above object as a compressed file so that it saves space and you don't have to regenerate it! Woot
save.lbzip2(l2.gex.qnorm, file = "rna/rObject_l2.gex.qnorm.RDataFS", n.cores = 16)

# With the above complete, subset the ATAC-seq object to only have the columns that the above RNA-seq file does
l2.atac.cpm.qnorm <- l2.atac.cpm.qnorm[, colnames(l2.atac.cpm.qnorm) %in% colnames(l2.gex.qnorm) ]; rownames(l2.atac.cpm.qnorm) <- gsub(".*\\|", "", rownames(l2.atac.cpm.qnorm))

# Get rid of the objects in memory you no longer need for space
rm(gex, l2.gex)

# To only operate on the correlations on the same chromosome, generate a list of each of the genes on each chromosome and put that into an object
# so you don't have to figure out which is on which chromosome every iteration...
i <- unique(seqnames(open)); pb <- lapply(i, function(j){subset(open, seqnames == j)$id}); names(pb) <- i
pc <- lapply(i, function(j){subset(open, seqnames == j)$coords}); names(pc) <- i

# Just generate all of the possible correlations you can across each chromosome to use as the expected correlations for each peak
d$coords <- with(d, paste0(V3, ":", start, "-", end))
genes <- rownames(l2.gex.qnorm)

# Get a list of the genes with multiple TSS sites for below
dup <- subset(as.data.frame(table(d$V13)), Freq > 1)$Var1

# First, generate the object with each combination of filtered peak and gene that are on the same chromosome! Yeah
cor.mat <- bind_rows( mclapply(1:length(genes), function(i){
	return( data.frame( "peak" = pb[[ subset(d, V13 == genes[i])[1,"V3"] ]], "gene" = genes[i], "peak.coords" = pc[[ subset(d, V13 == genes[i])[1,"V3"] ]], "gene.coords" = ifelse( genes[i] %in% dup, "multi", subset(d, V13 == genes[i])[1,"coords"] ) ) )
}, mc.cores = 16) )

# Remove the objects in memory you don't need anymore
rm(i, pb, pc)

# Use the power of mclapply to iterate through each of the above rows and calculate the distance between each peak and the gene TSS
a <- which( cor.mat$gene.coords == "multi" )
cor.mat[ a, "gene.coords" ] <- unlist( mclapply( a, function(i){ 
	ps <- as.numeric(gsub("-.*", "", gsub(".*:", "", cor.mat[i,"peak.coords"]))); pe <- as.numeric(gsub(".*-", "", cor.mat[i,"peak.coords"]))
	gs <- as.numeric(gsub("-.*", "", gsub(".*:", "", subset(d, V13 == cor.mat[i,"gene"])$coords))); ge <- as.numeric(gsub(".*-", "", subset(d, V13 == cor.mat[i,"gene"])$coords))
	b <- which( pe >= gs & ps <= ge )
	if ( length(b) != 0 ) { 
		return( d[J(cor.mat[i,"gene"]), nomatch = 0L]$coords[b[1]] ) 
	} else {
		tmp <- ifelse( pe <= gs & ps <= gs, pe - gs, ps - ge )
		return( d[J(cor.mat[i,"gene"]), nomatch = 0L]$coords[which(abs(tmp) == min( abs( tmp ) ))][1] )
	}
}, mc.cores = 64) )

# Turn the "cor.mat" object into a data.table to benefit from speedups granted to the object
cor.mat <- as.data.table(cor.mat)

cor.mat$dist <- unlist( mclapply(1:nrow(cor.mat), function(i){
	ps <- as.numeric(gsub("-.*", "", gsub(".*:", "", cor.mat[i,peak.coords]))); pe <- as.numeric(gsub(".*-", "", cor.mat[i,peak.coords]))
	gs <- as.numeric(gsub("-.*", "", gsub(".*:", "", cor.mat[i,gene.coords]))); ge <- as.numeric(gsub(".*-", "", cor.mat[i,gene.coords]))
	return( ifelse( pe >= gs & ps <= ge, 0, ifelse( pe <= gs & ps <= gs, pe - gs, ps - ge ) ) )
}, mc.cores = 32) ) 

# Remove all tested/generated objects in memory you no longer need
rm(a,d,dup,genes,open)

# Save the main object just in case something goes wrong
save.lbzip2(cor.mat, file = "atac/rObject_cor.mat.RDataFS", n.cores = 16)

# To attempt to speed things up, you changed the luts "l2.atac.cpm" and "l2.gex" into matrix objects
l2.atac.cpm.qnorm <- as.matrix(l2.atac.cpm.qnorm); l2.gex.qnorm <- as.matrix(l2.gex.qnorm)

# The object above has 461777254 rows (!) so break it up into 2 even pieces to minimize memory usage 
it <- data.frame("start" = ((0:1)*(nrow(cor.mat)/2))+1, "end" = (1:2)*(nrow(cor.mat)/2) )

# Save the above object to memory so that you can load it in and only take the rows you want
mclapply(1:2, function(i) {
	a <- cor.mat[c(it[i,"start"]:it[i,"end"]),]
	save.lbzip2(a, file = paste0("atac/rObject_cor_part_", i, ".mat.RDataFS"), n.cores = 30)
}, mc.cores = 2)


### NOTE: Here, and after each section, the environment was closed and reopened to minimize additional objects mucking up the analyses


## Section 1

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(fastSave); library(doParallel); library(data.table)
})

# Load in the current object you are working on. Note that the object will be named "a" 
load.lbzip2("atac/rObject_cor_part_1.mat.RDataFS", n.cores = 16)

# Start by loading in the log2 ATAC-seq data using CPM and then quantile normalized
load.lbzip2("atac/rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)
rownames(l2.atac.cpm.qnorm) <- gsub(".*\\|", "", rownames(l2.atac.cpm.qnorm))

# Now load in the quantile normalized log2 CPM RNA-seq data 
load.lbzip2("rna/rObject_l2.gex.qnorm.RDataFS", n.cores = 16)

# Filter the objects to minimize memory usage
l2.gex.qnorm <- as.matrix(subset(l2.gex.qnorm, rownames(l2.gex.qnorm) %in% a$gene))
l2.atac.cpm.qnorm <- as.matrix(subset(l2.atac.cpm.qnorm, rownames(l2.atac.cpm.qnorm) %in% a$peak)[ , colnames(l2.gex.qnorm) ])

# Time to generate the correlations
a$cor.s <- unlist( mclapply( 1:nrow(a), function(i) { return( suppressWarnings( cor( as.numeric(l2.atac.cpm.qnorm[a[i,peak],]), as.numeric(l2.gex.qnorm[a[i,gene],]), method = "spearman" ) ) ) }, mc.cores = 48 ) )

# Once it's done, save it to disk and remove all objects and get ready to start again
save.lbzip2( a, file = "atac/rObject_finished_cor_part_1.RDataFS", n.cores = 16 )


## Section 2

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(fastSave); library(doParallel); library(data.table)
})

# Load in the current object you are working on. Note that the object will be named "a" 
load.lbzip2("atac/rObject_cor_part_2.mat.RDataFS", n.cores = 16)

# Start by loading in the log2 ATAC-seq data using CPM and then quantile normalized
load.lbzip2("atac/rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)
rownames(l2.atac.cpm.qnorm) <- gsub(".*\\|", "", rownames(l2.atac.cpm.qnorm))

# Now load in the log2 RNA-seq data using CPM and then quantile normalized
load.lbzip2("rna/rObject_l2.gex.qnorm.RDataFS", n.cores = 16)

# Filter the objects to minimize memory usage
l2.gex.qnorm <- as.matrix(subset(l2.gex.qnorm, rownames(l2.gex.qnorm) %in% a$gene))
l2.atac.cpm.qnorm <- as.matrix(subset(l2.atac.cpm.qnorm, rownames(l2.atac.cpm.qnorm) %in% a$peak)[ , colnames(l2.gex.qnorm) ])

# Time to generate the correlations
a$cor.s <- unlist( mclapply( 1:nrow(a), function(i) { return( suppressWarnings( cor( as.numeric(l2.atac.cpm.qnorm[a[i,peak],]), as.numeric(l2.gex.qnorm[a[i,gene],]), method = "spearman" ) ) ) }, mc.cores = 48 ) )

# Once it's done, save it to disk and remove all objects and get ready to start again
save.lbzip2( a, file = "atac/rObject_finished_cor_part_2.RDataFS", n.cores = 16 )


## Bringing everything together

suppressPackageStartupMessages({
	library(fastSave); library(GenomicRanges); library(dplyr); library(data.table); library(doParallel)
})

# Load in each of the different objects into a single list object to combine them using bind_rows()
bl <- list()
load.lbzip2( "atac/rObject_finished_cor_part_1.RDataFS", n.cores = 16 )
bl[[1]] <- a
load.lbzip2( "atac/rObject_finished_cor_part_2.RDataFS", n.cores = 16 )
bl[[2]] <- a
bb <- as.data.table(bind_rows(bl)); rm(bl); gc()

# To take advantage of the power of data.table and its binary search for speedy accesses, set the key, which in this case are the peaks
setkey(bb, peak)

# Figure out what all of the thresholds are for each peak. Woot woot. Takes neither time nor memory due to the power of binary search
tmp <- unique(bb$peak)
thresh <- bind_rows(mclapply( tmp, function(i){ d <- bb[J(i), nomatch = 0L][!is.na(cor.s)]; return( data.frame( "top" = quantile(d$cor.s, .75) + (1.5*IQR(d$cor.s)), "bottom" = quantile(d$cor.s, .25) - (1.5*IQR(d$cor.s)), row.names = i ) ) }, mc.cores = 64))
bb$filt.s <- ifelse( ( thresh[bb$peak,"top"] <= bb$cor.s | thresh[bb$peak,"bottom"] >= bb$cor.s ) & abs(bb$dist) <= 1e6, T, F )

# Save the above object to not have to generate it again
write.table(thresh, file = "atac/rObject_thresh.txt", quote = F, sep = "\t")

# Start by adding the abc peak-gene linkages onto the bb object! But first, add the gene column to the keys so that it can be binary searched!
bb <- setkey(bb, peak, gene)
a <- arrange(read.table("microc/rss_abc_peak_int_all.txt.gz"), V11)
b <- unlist(mclapply(1:nrow(a), function(i){ bb[J(a[i,"V4"], a[i,"V9"]), which = T] }, mc.cores = 16))
bb$abc <- 0.0; bb$filt.abc <- F
bb[b,"abc"] <- a$V11; bb[b,"filt.abc"] <- T

# Load in the list of peak/accessible te to split the correlation measures and attempt to see what's up with ODE elements overall
int <- read.table("atac/peak.te.int"); pte <- unique(int$V4); pnonte <- tmp[!tmp %in% pte]

# Annotate the peaks as being TE-derived, nonTE-derived, or ODE-derived!
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
pte <- pte[!pte %in% ode$peak]; r1e <- subset(ode, repName == "ORR1E")$peak; r1d2 <- subset(ode, repName == "ORR1D2")$peak; r1d2 <- r1d2[!r1d2 %in% r1e]

# Annotate the peaks as ORR1E or ORR1D2. In the cases of a peak containing both, annotate it as ORR1E
bb$g1 <- ""
for (f in pte) { set(bb, bb[J(f), which = T] ,"g1","te") }
for (f in pnonte) { set(bb, bb[J(f), which = T] ,"g1","nonte") }
for (f in r1e) { set(bb, bb[J(f), which = T] ,"g1","orr1e") }
for (f in r1d2) { set(bb, bb[J(f), which = T] ,"g1","orr1d2") }

# Now annotate the clusters for ORR1E and ORR1D2
bb$g2 <- bb$g1
for (f in r1e) { set(bb, bb[J(f), which = T], "g2", subset(ode, peak == f & repName == "ORR1E")$cpm.c[1] ) }
for (f in r1d2) { set(bb, bb[J(f), which = T], "g2", subset(ode, peak == f & repName == "ORR1D2")$cpm.c[1] ) }

# Save the object to not have to regenerate it again, overwritting the previous one in the process
cor.mat <- bb
save.lbzip2(cor.mat, file = "atac/rObject_cor.mat.RDataFS", n.cores = 16)

# Remove all the objects you no longer need
rm(a,b,bb,f,int,pnonte,pte,r1d2,r1e,thresh,tmp)






