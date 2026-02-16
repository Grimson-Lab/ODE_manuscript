# This script will iterate through the output of blast and attempt to merge nearby blast hits for the same element that might be split by indel(s) of some kind. 
# If the blast row has most, ~85-90%, of its sequence found, skip and go to the next one. Otherwise, keep looking for nearby chunks to combine together
# 
# Expected script call (post running module load R/4.1.3-r9):
# Rscript --vanilla blastMerge.r [input.blast] [output.blast] 

# Load in the necessary libraries
suppressPackageStartupMessages({ require(doParallel); library(dplyr) }) 

# Load in the args that were used to initiate the script 
args = commandArgs(trailingOnly = T)

# Load in the input file - initial output from blast for a given strain
ib <- read.table(args[1])

# Set up the list of element IDs to iterate through and then do just that! Append the values to a print data.frame to save as a file later
samps <- unique(ib$V1); out <- data.frame(matrix(nrow = length(samps), ncol = ncol(ib)))

# Time for the main method! Using mclapply makes things go brrrrrr
results <- mclapply(1:length(samps), function(f) {
	sub <- subset(ib, V1 == samps[f])

	# If there is only one entry in sub, ez pz. Otherwise, iterate over each of the rows and attempt to find the "best" merged interval. Use the combined 
	# bitscore produced by blast for finding the "best" one
	if ( nrow(sub) == 1 ) { out[f,] <- sub[1,]; return(out[f,]); break }

	# You wanted to make a variable to discourage merges, which is applied when combining bitscores (column 13)
	penalty <- .95

	# Initialize the variables for the "current" entry and set the "best" result to the first row of sub 
	count <- 1; cC <- sub[1,2]; cS <- sub[1,9]; cE <- sub[1,10]; cQS <- sub[1,7]; cQE <- sub[1,8]; cB <- sub[1,13]; cSd <- sub[1,14]; m <- 1
	out[f,] <- sub[1,]; tmp <- out[f,]
	for (g in 1:(nrow(sub)-1) ) {

		# Initialize the variables for the "next" entry
		nC <- sub[g+1,2]; nS <- sub[g+1,9]; nE <- sub[g+1,10]; nQS <- sub[g+1,7]; nQE <- sub[g+1,8]; nB <- sub[g+1,13]; nSd <- sub[g+1,14]

		# If m is 1, check whether the next interval is + or - relative to the current interval. Otherwise, check to make sure the strands are the same,
		# the intervals are close enough together to be merged, and that there isn't a large fraction of interval overlap
		if ( m == 1 ) { 
			if ( nC == cC & abs((nS-cE)-(nQS-cQE)) <= ((cE-cS+nE-nS)*.15) & (cE-nS < (nE-nS)*.3) & cSd == nSd ) { pass <- T } else { pass <- F }
		} else { 
			if ( nC == cC & abs((nS-cE)-(nQS-cQE)) <= ((cE-cS+nE-nS)*.15) & (cE-nS < (nE-nS)*.3) & tmp[1,14] == nSd ) { pass <- T } else { pass <- F } 
		}

		# Check to make sure the chromosomes match, the intervals are not too far apart, and their matched sequences could be reasonable to merge. If 
		# this fails, compare against the current best and swap if its bitscore is higher. Make sure to have a case for the final row if it fails to merge
		if ( pass ) {

			# Have the appended values be stored in this tmp object that will get added to for multiple instances as needed
			if ( m == 1 ) {
				tmp[1,] <- c(sub[g,1:3], sub[g,4:6] + sub[g+1,4:6], sub[g,7], sub[g+1,8], sub[g,9], sub[g+1,10:11], min(sub[g,12], sub[g+1,12]), (sub[g,13] + sub[g+1,13]) * penalty, sub[g+1,14])
			} else {
				tmp[1,] <- c(tmp[1,1:3], tmp[1,4:6] + sub[g+1,4:6], tmp[1,7], sub[g+1,8], tmp[1,9], sub[g+1,10:11], min(tmp[1,12], sub[g+1,12]), tmp[1,13] + (sub[g+1,13] * penalty), tmp[1,14])
			}

			# Increment m so that it is no longer zero and future merges will use tmp rather than sub
			m <- m + 1
		} else { 
			# If m is 1, use sub to check for a higher match. Otherwise, use tmp
			if ( m == 1 ) {
				# Check to see if the current bitscore is more than the current "best" match. If yes, replace
				if ( sub[g,13] > out[f,13] ) { out[f,] <- sub[g,] }
			} else {
				# Check to see if the current bitscore is more than the current "best" match. If yes, replace
				if ( (tmp[1,13] * penalty) > out[f,13] ) { out[f,] <- tmp[1,] }
			}

			# Case for the final row failing to merge and potentially being the "best" match 
			if ( g == (nrow(sub)-1) ) { if ( sub[g+1,13] > out[f,13] ) { out[f,] <- sub[g+1,] } }

			# Set m back to 1 and reset tmp (?) since a row merge failed 
			m <- 1
		}

		# Setting the current variables to their next values
		cC <- sub[g+1,2]; cS <- sub[g+1,9]; cE <- sub[g+1,10]; cQS <- sub[g+1,7]; cQE <- sub[g+1,8]; cB <- sub[g+1,13]; cSd <- sub[g+1,14]

		# Case for if the final row merges and potentially being the "best" match
		if ( g == (nrow(sub)-1) ) { if ( tmp[1,13] > out[f,13] ) { out[f,] <- tmp[1,] } }
	}

	# Now return the out[f,] object so that it can be combined after mclapply for speedup purposes
	return(out[f,])

}, mc.cores = 64)

# Output the file so it can be used with the awk filter command into bedtools intersect
write.table(bind_rows(results), args[2], sep = "\t", quote = F, row.names = F, col.names = F)
