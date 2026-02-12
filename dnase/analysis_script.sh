# !/bin/bash

# This script will download all of the DNase-seq narrowPeak files from the Mouse Developmental Matrix project from ENCODE so that you can run
# your TEshuffle script and get a sense of whether ODEs are significantly enriched in bulk-tissue accessibility. Hurrah!  

######################################
### Download the data!
######################################

# State the current task
echo "#######################"
echo "Downloading raw data..."
echo "#######################"

# Using the ENCODE portal, download all the "narrowPeak" files for DNase-seq experiments from the Mouse Developmental Matrix 
wget "https://www.encodeproject.org/metadata/?status=released&related_series.%40type=OrganismDevelopmentSeries&replicates.library.biosample.organism.scientific_name=Mus+musculus&assay_title=DNase-seq&files.file_type=bed+narrowPeak&type=Experiment&files.analyses.status=released&files.preferred_default=true" -O dnase/DNase.metadata

# Get a list of all the commands you'll need to run to download all the files, and then download them using "parallel"
# mkdir dnase/0_beds
awk -F'\t' '{gsub(" ","_",$11); if ($2 == "bed narrowPeak") print "wget -q " $48 " -O dnase/0_beds/" $11}' dnase/DNase.metadata | awk -F' ' '{x[$5]++; print $0 "_" x[$5] ".bed.gz"}' > dnase/wget_cmds
parallel -j 24 < dnase/wget_cmds

# State the current task
echo "################################"
echo "Finished downloading raw data..."
echo "################################"



######################################
### Shuffle for TE enrichment! 
######################################

# State the current task
echo "######################"
echo "Starting TE Shuffle..."
echo "######################"

# Use the "TEshuffle.r" script to take the peaks not overlapping coding sequence and measure enrichment of TEs across each of the above datasets
module load R/4.1.3-r9
Rscript scripts/TEshuffle.r dnase/0_beds dnase/1_shuffleResults *.bed.gz 40

# Now that the TE Shuffle is done and has summarized all the results, get rid of each individual file
rm dnase/1_shuffleResults/.*

# State the current task
echo "######################"
echo "Finished TE Shuffle..."
echo "######################"


