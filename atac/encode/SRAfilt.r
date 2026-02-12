# This script will load in a given file downloaded from the SRA Run Selector and turn it into a tab-deliniated file for easier analyses
# 
# Expected script call:
# Rscript --vanilla SRAfilt.r [SRA Run Selector file] [output file name]

# Load in the file name arguement
args <- commandArgs(trailingOnly = T)

# Read in the file
a <- subset(read.csv(args[1]), dev_stage %in% c("embryonic", "postnatal") & sample_type == "tissue")

# Below includes manual changes based on the experiment and shortening particular information (such as assay, antibody name, etc)
a$id <- with(a, ifelse( dev_stage == "embryonic", paste0( gsub(" ", "_", tissue), "_e", gsub(" day", "", AGE) ), paste0( gsub(" ", "_", tissue), "_p", gsub(" day", "", AGE) ) ) )

# Finally, output the file and only keep the columns you care about
write.table( a[ , c("Run", "AvgSpotLen", "id") ], file = args[2], row.names = F, sep = "\t", quote = F )
