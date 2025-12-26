# !/bin/bash

# This script encapsulates all of your RNA-seq pipelines for paired-end data using a manually added SRA Run Selector list of data/labels. 
# Use the first and second parameters to specify the specific "SRAfilt" R script, the input to that script, and the output directory for 
# everything!
# 
# Expected script call
# sh analysis_script [output directory]


########################
### Run RepeatMasker ###
########################

# State the current task
echo "########################"
echo "Starting RepeatMasker..."
echo "########################"

# Run RepeatMasker on mm10 (will take a few days...)
~/RepeatMasker/RepeatMasker -species "mus musculus" -s -e hmmer -pa 30 -no_is -nolow -dir $1 mm10.fa

# State the current task
echo "########################"
echo "Finished RepeatMasker..."
echo "########################"



# State the current task
echo "################################"
echo "Starting merging and trimming..."
echo "################################"

# Run the python script to merge nearby TE annotations that should be labeled as a single element
python "${1}"/controlMerge.py "${1}"/Genome.out "${1}"/merge.txt

# Run the python script to remove all TE annotations <50bps and prevent overlapping annotations
python "${1}"/coordTrim.py "${1}"/merge.txt 50 "${1}"/mm10_repeatMasker.txt

# State the current task
echo "################################"
echo "Finished merging and trimming..."
echo "################################"


# Compress the above output and remove all unneccesary files
pigz -p 16 "${1}"/mm10_repeatMasker.txt
rm "${1}"/Genome.* "${1}"/merge.txt


