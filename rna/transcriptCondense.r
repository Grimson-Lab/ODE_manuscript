# This script will read in the output of StringTie across all ImmGen RNA-seq datasets as well as the human-readable names and condense each set
# of transcripts into a single row for future analyses
# 
# Expected script call:
# Rscript --vanilla transcriptCondense.r [StringTie output directory] [pattern to match] [output file name]

# Load in the necessary libraries 
suppressPackageStartupMessages({
	library(doParallel); library(dplyr)
})

# Load in the file name arguement
args <- commandArgs(trailingOnly = T)

# If the path argument doesn't end in a "/", make it so
if ( !grepl( "/$", args[1] ) ) { args[1] <- paste0(args[1], "/") }

# Make sure to change these god awful names to their respective symbol. To do this, load in the file with that info!
ref <- subset(read.table("../ncbiRefSeqCurated.txt.gz", sep = "\t"), !grepl("_alt|chrUn|_random|_fix|chrY|chrM", V3) )
rownames(ref) <- ref$V2

# Identify all of the samples you'll need to load in and then read in each of the respective files, condensing the TPM values of all 
# the different isoforms into a single row with the gene symbol rather than the Refseq ID
n <- list.files(path = args[1], pattern = args[2])
f <- bind_cols( mclapply(n, function(i){
	a <- read.table( paste0(args[1], i), header = T, sep = "\t" ); colnames(a) <- c("rsID", "gn", "chr", "strand", "start", "end", "cov", "fpkm", "tpm")
	a$gn <- ref[a$rsID, "V13"] 
	a <- a %>% group_by(gn) %>% summarize(tmp = sum(tpm)) %>% as.data.frame() %>% arrange(gn); colnames(a)[2] <- gsub(".out", "", i)
	rownames(a) <- a$gn

	return( a[ , gsub(".out", "", i), drop = F ] )
}, mc.cores = 16) )

### NOTE: Since stringtie is using the same gene annotations for each sample, the rownames will be identical  

# Finally, output the finished file
write.table(f, file = args[3], sep = ",", quote = F)


