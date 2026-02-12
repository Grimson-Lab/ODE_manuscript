# This script will load in the manually edited ImmGen RNA-seq metadata file downloaded from the SRA Run Selector 
# (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA429735&o=acc_s%3Aa) 
# and output a shortened file used for downloading and analyzing the data, specifically the samples that also exist in the ATAC-seq
# 
# Expected script call:
# Rscript --vanilla SRAfilt_encode.r [SRA Run Selector file] [output file name]

# Load in the file name arguement
args <- commandArgs(trailingOnly = T)

# Read in the file
rna <- read.table(args[1], header = T, sep = "\t")

# Read in the ATAC-seq metadata file
atac <- read.table("atac/immgen/SraRunTable_immgenATAC.tab", header = T, sep = "\t")

# Identify which RNA-seq samples are present in the ATAC-seq
b <- with(rna, which( gsub("_[1-9]", "", label) %in% gsub("_[1-9]", "", atac$label) ))

# Finally, output the file with each of the samples you want and only keep the columns you care about
write.table(rna[ b, c("Run", "AvgSpotLen", "label") ], file = args[2], sep = "\t", quote = F, row.names = F)

