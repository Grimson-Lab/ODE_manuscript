### Tuesday January 28th, 2025
# This script will download all of the DNase-seq narrowPeak files from the Mouse Developmental Matrix project from ENCODE so that you can run
# your TEshuffle script and get a sense of whether ODEs are significantly enriched in bulk-tissue accessibility. Hurrah!  

######################################
### Download the data!
######################################

# Using the ENCODE portal, download all the "narrowPeak" files for DNase-seq experiments from the Mouse Developmental Matrix 
wget "https://www.encodeproject.org/metadata/?status=released&related_series.%40type=OrganismDevelopmentSeries&replicates.library.biosample.organism.scientific_name=Mus+musculus&assay_title=DNase-seq&files.file_type=bed+narrowPeak&type=Experiment&files.analyses.status=released&files.preferred_default=true" -O DNase.metadata

# Get a list of all the commands you'll need to run to download all the files, and then download them using "parallel"
mkdir 0_beds
awk -F'\t' '{gsub(" ","_",$11); if ($2 == "bed narrowPeak") print "wget -q " $48 " -O 0_beds/" $11}' DNase.metadata | awk -F' ' '{x[$5]++; print $0 "_" x[$5] ".bed.gz"}' > wget_cmds
parallel -j 24 < wget_cmds


######################################
### Shuffle for TE enrichment! 
######################################

# Use the "TEshuffle.r" script to take the peaks not overlapping coding sequence and measure enrichment of TEs across each of the above datasets
module load R/4.1.3-r9
Rscript /workdir/jdc397/1_currentWork/8_teCactus/00_clean/1_data/1_atac/TEshuffle.r 0_beds 1_shuffleResults *.bed.gz 40




