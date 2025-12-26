# This script will load in the manually edited ImmGen ATAC-seq metadata file downloaded from the SRA Run Selector and output a shortened file 
# used for downloading and analyzing the data
# 
# Expected script call (post running module load R/4.1.3-r9):
# Rscript --vanilla SRAfilt.r [SRA Run Selector file] [output file name]

# Load in the file name arguement
args <- commandArgs(trailingOnly = T)

# Read in the file
a <- read.table(args[1], header = T, sep = "\t")

# Finally, output the file and only keep the columns you care about
write.table(a[ a$AvgSpotLen == "75", c("Run", "AvgSpotLen", "label")], file = args[2], sep = "\t", quote = F, row.names = F)
