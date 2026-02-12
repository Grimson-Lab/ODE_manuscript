# !/bin/bash

# This script encapsulates all of your RNA-seq pipelines for paired-end data using a manually added SRA Run Selector list of data/labels. 
# Use the first and second parameters to specify the specific "SRAfilt" R script, the input to that script, and the output directory for 
# everything!
# 
# Expected script call (from "ODE_manuscript directory")
# sh analysis_script [R code name] [R code input] [R code output] [output directory] 

### NOTE: For the above data, you downloaded the metadata from the SRA Run Selector (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA429735&o=acc_s%3Aa)
### and then *manually* linked the SRA numbers to the sample names from ImmGen in their Yoshida et al 2019 paper. That file is "SraRunTable_immgenRNA.tab"

#########################################
### Download the Data (with Parallel!)
#########################################

# Start by making the directory to put the files in *if* it doesn't already exist
mkdir -p $4

# Use an R script to load in the SRA run selector datasets to then use awk to create all the run files
Rscript --vanilla "${4}"/"${1}" "${4}"/"${2}" "${4}"/"${3}"

# Automate picking a number of cores to use that won't crash the server... "cores" is the number of parallel jobs, and "cores2" is the number
# of cores per each individual job
cores=`sed '1d' "${4}"/"${3}" | wc -l`; cores2=8
if [ $cores -gt 6 ]
then 
	cores=6
fi

# State the current task
echo "#######################"
echo "Downloading raw data..."
echo "#######################"

# Now that you have your file, generate the commands to get all of the files downloaded using awk and fasterq-dump for the paired-end data
awk -v a=$cores2 -v b=$4 '{FS="\t"; OFS=""; if ( NR != 1 ) { print "fasterq-dump ", $1 " -S; pigz -c -p " a " " $1 "_1.fastq > " b "/" $3 "a.fq.gz; pigz -c -p " a " " $1 "_2.fastq > " b "/" $3 "b.fq.gz; rm " $1 "_1.fastq " $1 "_2.fastq" } }' < "${4}"/"${3}" > "${4}"/dump_cmds
parallel -j $cores < "${4}"/dump_cmds

# State the current task
echo "################################"
echo "Finished downloading raw data..."
echo "################################"



#####################################
### Trim the Data (with Parallel!)
#####################################

# Make the output directories
mkdir "${4}"/trim; mkdir "${4}"/trim/fastqc

# Make sure to add the most recent cutadapt path which uses python3.0+ and allows for multi-threading, which in our server is the below:
export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages
export PATH=/programs/cutadapt-4.1/bin:$PATH

# State the current task
echo "########################"
echo "Starting the trimming..."
echo "########################"

# Write the trim commands to a file and then run those commands using parallel. Check to make sure the length is approximately 25% of the total
# read size (usually given as full read size, not the size of each mate), using the power of variables and math
awk -v a=$cores2 -v b=$4 '{FS="\t"; OFS=""; if ( NR != 1 ) { print "trim_galore -j " a " --stringency 4 --length " int($2 / 8 + 0.5) " --paired -o " b "/trim " b "/" $3 "a.fq.gz " b "/" $3 "b.fq.gz --fastqc_args \"--outdir " b "/trim/fastqc/\" 2>/dev/null"} }' < "${4}"/"${3}" > "${4}"/trim_cmds
parallel -j $cores < "${4}"/trim_cmds

# State the current task
echo "####################"
echo "Finished trimming..."
echo "####################"



######################################
### Align the Data (with Parallel!)
######################################

# Load the path needed for one of the recent versions of STAR so that you can randomly choose a read to be used for multi-mapping purposes
export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH

# Create the directory where the alignments will live 
mkdir "${4}"/STAR

# State the current task
echo "#####################"
echo "Starting alignment..."
echo "#####################"

# Check to make sure the .gtf file is uncompressed for the below STAR command
if [ -f mm10.refseqCurated.gtf.gz ]; then gunzip mm10.refseqCurated.gtf.gz; fi

# Populate the command file and then align using parallel 
> "${4}"/star_cmds
for f in $(ls "${4}"/trim/*val*.gz | sed 's/[a-b]_val_[1-2].fq.gz//' | uniq)
do 
	n=`echo $f | sed 's|.*/||g'`
	### NOTE: The below will throw a syntax error but still works to do the correct math. *shrug*
	l=$(grep `echo $n` "${4}"/SRA.tab | cut -f2); m=$(expr $l / 2 - 1)
	echo "STAR --genomeDir STAR_noGene_mm10 --outFileNamePrefix "${4}"/STAR/"${n}" --sjdbGTFfile mm10.refseqCurated.gtf --sjdbOverhang "${m}" --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFilterMatchNminOverLread .4 --seedSearchStartLmax 15 --readFilesCommand zcat --outFilterMultimapNmax 500 --outSAMmultNmax 1 --outMultimapperOrder Random --readFilesIn "${f}"a_val_1.fq.gz "${f}"b_val_2.fq.gz --runMode alignReads --alignEndsType Local --runThreadN "${cores2}"" >> "${4}"/star_cmds
done
parallel -j $cores < "${4}"/star_cmds

# State the current task
echo "#####################"
echo "Finished alignment..."
echo "#####################"


# State the current task
echo "##############################"
echo "Starting merging & deduping..."
echo "##############################"

# Check to see if you have any replicates for the paired-end reads. If not, don't merge. Otherwise, merge

# Now that you have the aligned data, merge all of the replicates (if present) and then run Picard to remove the PCR/Optical duplicates as 
# well as filter out all of the reads that align to the gross blacklist regions and chrM/Y
> "${4}"/picard_cmds
for f in $(ls "${4}"/STAR/*Aligned.sortedByCoord.out.bam | sed 's/_[1-9]Aligned.*//g' | uniq)
do 
	if [[ `ls $f*Aligned* | wc -l` -gt 1 ]]
	then 
		echo "samtools merge -@ "${cores2}" -o "$f"_merge.bam $(ls "${f}"*.out.bam | tr '\n' ' '); samtools view -@"${cores2}" -hf 2 -q 10 "$f"_merge.bam | awk '{if (\$3 != \"chrM\" && \$3 != \"chrY\") print}' | samtools view -b -@"${cores2}" | bedtools intersect -v -wa -a stdin -b comb_badRegions.bed > "$f"_tmp.bam; java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates I="$f"_tmp.bam ASSUME_SORT_ORDER=coordinate O="$f".dedup.bam M="$f".metrics REMOVE_DUPLICATES=true" >> "${4}"/picard_cmds
	else 
		# For any data with just one replicate, thus no merge needed
		echo "samtools view -@"${cores2}" -hf 2 -q 10 $(ls "${f}"*.out.bam) | awk '{if (\$3 != \"chrM\" && \$3 != \"chrY\") print}' | samtools view -b -@"${cores2}" | bedtools intersect -v -wa -a stdin -b comb_badRegions.bed > "$f"_tmp.bam; java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates I="$f"_tmp.bam ASSUME_SORT_ORDER=coordinate O="$f".dedup.bam M="$f".metrics REMOVE_DUPLICATES=true" >> "${4}"/picard_cmds
	fi
done
parallel -j $cores < "${4}"/picard_cmds

# State the current task
echo "##############################"
echo "Finished merging & deduping..."
echo "##############################"

# Now remove all the original .fq.gz files and trimmed .fq.gz files
rm "${4}"/*.fq.gz "${4}"/trim/*.fq.gz "${4}"/STAR/*.tab  



######################################
### Stringtie the Data (with Parallel!)
######################################

# Make the output directory
mkdir "${4}"/stringtie

# Load the path needed for stringtie
export PATH=/programs/stringtie-2.2.1:$PATH

# State the current task
echo "####################################"
echo "Starting Stringtie quantification..."
echo "####################################"

# Use this to write the commands to a file and then run those commands using parallel with 8 jobs with 3 cores each?
> "${4}"/stringtie_cmds
for f in "${4}"/STAR/*.bam
do 
	n=`echo $f | sed 's|.*/||g' | sed 's/.dedup.bam//g'`
	echo "stringtie -p 8 -f 0.1 -m 200 -a 10 -j 1 -g 50 -M 0.95 -c 2.5 -e -A "${4}"/stringtie/"${n}".out -o "${4}"/stringtie/"${n}".gtf -G mm10.refseqCurated.gtf $f 2>/dev/null" >> "${4}"/stringtie_cmds
done
parallel -j $cores < "${4}"/stringtie_cmds

# Combine all the TPM values across each sample together using the below R script and change the 
# names to be something more meaningful than Refseq IDs
Rscript --vanilla transcriptCondense.r "${4}"/stringtie/ .*out refseqCurated_tpm.csv

# With the gene counts done, remove everything within the trim, STAR, and stringtie folders you won't need anymore
rm -R "${4}"/trim "${4}"/STAR "${4}"/stringtie

# State the current task
echo "####################################"
echo "Finished Stringtie quantification..."
echo "####################################"





