## This file contains all the necessary code to download human RNA-seq data in various immune cells to use as a comparison to the Immgen mouse
## data 

# Make sure to add "lbzip2" to your environment prior to starting R
# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH

# Load in the packages you'll need throughout this section
suppressPackageStartupMessages({ 
	library(fastSave); library(data.table); setDTthreads(8); library(dplyr)
}) 

# Start with downloading data in 18 cell types *and* total PBMCs
system("wget https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip")

# Get the data for 30 immune cell types from Monaco et al
system("wget https://www.proteinatlas.org/download/tsv/rna_immune_cell_monaco.tsv.zip")

# Finally, download the data for 15 immune cell types from Schmiedel BJ et al 2018 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6289654/)
system("wget https://www.proteinatlas.org/download/tsv/rna_immune_cell_schmiedel.tsv.zip")

# Uncompress them all! 
system("for f in *.zip; do unzip $f; done")

# Time to load in each of the files! 
s <- as.data.frame(fread("rna_immune_cell_schmiedel.tsv")); s$data <- "schmiedel"
m <- as.data.frame(fread("rna_immune_cell_monaco.tsv")); m$data <- "monaco"
b <- as.data.frame(fread("rna_immune_cell.tsv")); b$data <- "hpa"

# You want to take cell types that are representative of the datasets that you have from ImmGen
x <- rbind(s,m[,c(1:4,6)],b[,c(1:4,7)]); colnames(x) <- c("id", "name", "type", "tpm", "data")

a <- c("NK-cell" = "NK", "classical monocyte" = "C.Mo", "naive CD8 T-cell" = "CD8", "myeloid DC" = "mDC", "naive B-cell" = "B", "naive CD4 T-cell" = "CD4", 
	"Naive CD4 T-cell activated" = "Act.CD4", "Naive CD8 T-cell activated" = "Act.CD8", "intermediate monocyte" = "I.Mo", "non-classical monocyte" = "NC.Mo", 
	"Vd2 gdTCR" = "gdT", "plasmacytoid DC" = "pDC", "gdT-cell" = "gdT", "memory B-cell" = "Mem.B", "memory CD4 T-cell" = "Mem.CD4", "memory CD8 T-cell" = "Mem.CD8", 
	"basophil" = "basophil", "eosinophil" = "eosinophil")
x <- subset(x, type %in% names(a)); x$short.type <- a[x$type]; x$comb <- paste0(x$id, "|", x$name); x$comb2 <- paste0(x$short.type, "_", x$data)
x <- arrange(x, short.type)

# Create the underlying matrix that you'll use, and populate it with all the gene information!
xm <- as.data.frame(matrix(0, nrow = length(unique(x$comb)), ncol = length(unique(x$comb2))), row.names = unique(x$comb)); colnames(xm) <- unique(x$comb2)
for (f in 1:nrow(x)) {
	xm[ x[f,"comb"], x[f,"comb2"] ] <- x[f,"tpm"]
}

# Output the new file version
write.table(log2(xm+1), file = "rObject_l2.hpa.txt", sep = "\t", quote = F)

# Perform qnorm on the above object
system("qnorm rObject_l2.hpa.txt > rObject_l2.hpa.qnorm.txt")

# Load the qnorm output 
l2.hpa.qnorm <- read.table("rObject_l2.hpa.qnorm.txt", header = T, row.names = 1)

# Save the above object as a compressed file so that it saves space and you don't have to regenerate it! Woot
save.lbzip2(l2.hpa.qnorm, file = "rObject_l2.hpa.qnorm.RDataFS", n.cores = 16)

# Remove all the files you no longer need
file.remove("rObject_l2.hpa.qnorm.txt")
system("rm *tsv*")
