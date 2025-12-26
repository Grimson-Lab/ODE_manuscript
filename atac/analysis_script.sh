# !/bin/bash

# This script encapsulates all of your ATAC-seq pipelines for paired-end data using a manually added SRA Run Selector list of data/labels. 
# Use the first and second parameters to specify the specific "SRAfilt" R script, the input to that script, and the output directory for 
# everything!
# 
# Expected script call
# sh analysis_script [R code name] [R code input] [R code output] [output directory] 

######################################
### Download the Data
######################################

# Start by making the directory to put the files *if* it doesn't already exist
mkdir -p $4

# Use an R script to load in the SRA run selector datasets to then use awk to create all the run files
Rscript --vanilla $1 $2 "${4}"/"${3}" 

# Automate picking a number of cores to use that won't crash the server... "cores" is the number of parallel jobs, and "cores2" is the number
# of cores per each individual job
cores=`sed '1d' "${4}"/"${3}" | wc -l`; cores2=8
if [ $cores -gt 6 ]
then 
	cores=6
fi

# Now that you have your file, generate the commands to get all of the files downloaded using awk and fasterq-dump for the paired-end data
awk -v a=$cores2 -v b=$4 '{FS="\t"; OFS=""; c[$3]++; if ( NR != 1 ) { print "fasterq-dump ", $1 " -S; pigz -c -p " a " " $1 "_1.fastq > " b "/" $3 "_" c[$3] "a.fq.gz; pigz -c -p " a " " $1 "_2.fastq > " b "/" $3 "_" c[$3] "b.fq.gz; rm " $1 "_1.fastq " $1 "_2.fastq" } }' < "${4}"/"${3}" > "${4}"/dump_cmds
parallel -j $cores < "${4}"/dump_cmds


######################################
### Trim the Data (with Parallel!)
######################################

# Make the output directories
mkdir "${4}"/trim; mkdir "${4}"/trim/fastqc

# Make sure to add the most recent cutadapt path which uses python3.0+ and allows for multi-threading!
export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages
export PATH=/programs/cutadapt-4.1/bin:$PATH

# Write the trim commands to a file and then run those commands using parallel. Check to make sure the length is approximately 25% of the total
# read size, depending on vibes. Luckily you can automate this now, so yay
awk -v a=$cores2 -v b=$4 '{FS="\t"; OFS=""; c[$3]++; if ( NR != 1 ) { print "trim_galore -j " a " --stringency 4 --length " int($2 / 8 + 0.5) " --paired -o " b "/trim " b "/" $3 "_" c[$3] "a.fq.gz " b "/" $3 "_" c[$3] "b.fq.gz --fastqc_args \"--outdir " b "/trim/fastqc/\" 2>/dev/null"} }' < "${4}"/"${3}" > "${4}"/trim_cmds
parallel -j $cores < "${4}"/trim_cmds


######################################
### Align the Data (with Parallel!)
######################################

# Load the path needed for one of the recent versions of STAR so that you can randomly choose a read to be used for multi-mapping purposes
export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH

# Create the directory where the alignments will live 
mkdir "${4}"/STAR

# Populate the command file and then align using parallel 
> "${4}"/star_cmds
for f in $(ls "${4}"/trim/*val*.gz | sed 's/[a-b]_val_[1-2].fq.gz//' | uniq)
do 
	n=`echo $f | sed 's|.*/||g'`
	echo "STAR --genomeDir ../STAR_noGene_mm10.vM25 --outFileNamePrefix "${4}"/STAR/"${n}" --alignMatesGapMax 2000 --outSAMtype BAM SortedByCoordinate --outFilterMatchNminOverLread .4 --seedSearchStartLmax 15 --readFilesCommand zcat --outFilterMultimapNmax 500 --outSAMmultNmax 1 --outMultimapperOrder Random --alignIntronMax 1 --readFilesIn "${f}"a_val_1.fq.gz "${f}"b_val_2.fq.gz --runMode alignReads --alignEndsType Local --runThreadN "${cores2}"" >> "${4}"/star_cmds
done
parallel -j $cores < "${4}"/star_cmds

# Check to see if you have any replicates for the paired-end reads. If not, don't merge. Otherwise, merge

# Now that you have the aligned data, merge all of the replicates and then run Picard to remove the PCR/Optical duplicates as well as filter
# out all of the gross blacklist reads
> "${4}"/picard_cmds
for f in $(ls "${4}"/STAR/*Aligned.sortedByCoord.out.bam | sed 's/_[1-9]Aligned.*//g' | uniq)
do 
	if [[ `ls $f*Aligned* | wc -l` -gt 1 ]]
	then 
		echo "samtools merge -@ "${cores2}" -o "$f"_merge.bam $(ls "${f}"*.out.bam | tr '\n' ' '); samtools view -@"${cores2}" -hf 2 -q 10 "$f"_merge.bam | awk '{if (\$3 != \"chrM\" && \$3 != \"chrY\") print}' | samtools view -b -@"${cores2}" | bedtools intersect -v -wa -a stdin -b ../comb_badRegions.bed > "$f"_tmp.bam; java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates I="$f"_tmp.bam ASSUME_SORT_ORDER=coordinate O="$f".dedup.bam M="$f".metrics REMOVE_DUPLICATES=true" >> "${4}"/picard_cmds
	else 
		# Just for the Stem cell data since it only has one replicate each
		echo "samtools view -@"${cores2}" -hf 2 -q 10 $(ls "${f}"*.out.bam) | awk '{if (\$3 != \"chrM\" && \$3 != \"chrY\") print}' | samtools view -b -@"${cores2}" | bedtools intersect -v -wa -a stdin -b ../comb_badRegions.bed > "$f"_tmp.bam; java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates I="$f"_tmp.bam ASSUME_SORT_ORDER=coordinate O="$f".dedup.bam M="$f".metrics REMOVE_DUPLICATES=true" >> "${4}"/picard_cmds
	fi
done
parallel -j $cores < "${4}"/picard_cmds

# Remove the temporary files for paired-end data
rm "${4}"/STAR/*_tmp.bam "${4}"/STAR/*Aligned.sortedByCoord.out.bam "${4}"/STAR/*_merge.bam

# Now remove all the original .fq.gz files and the intermediate .bam files you no longer need
rm "${4}"/*.fq.gz; rm "${4}"/STAR/*.tab


######################################
### Macs3 the Data 
######################################

# Make the output directory
mkdir "${4}"/macs3 

# Turn the paired-end file into individual cut sites and shift them based on the expected overhangs from Tn5 cutting
> "${4}"/bed_cmds
for f in $(ls "${4}"/STAR/*dedup.bam)
do
	echo "bamToBed -i $f | awk '{OFS=\"\\t\"; if (\$6 == \"+\") {print \$1, \$2 + 5, \$3, \$4, \$5, \$6} else {print \$1, \$2, \$3 - 5, \$4, \$5, \$6} }' > "${f:0:-9}"bed" >> "${4}"/bed_cmds
done
parallel -j $cores < "${4}"/bed_cmds

# Use the power of MACS3 to call peaks using parameters that follow the recommended ATAC-seq analysis pipeline 
source ~/Macs3Env/bin/activate 

# Time to run MACS3. You decided against using the summits since there were so many peaks otherwise and they were fairly small sizes
> "${4}"/macs3_cmds
for f in $(ls "${4}"/STAR/*.bed)
do 
	n=`echo $f | sed 's|.*/||g'`
	echo "macs3 callpeak -t $f -f BED -g mm -n "${n:0:-4}" -q 0.001 --nomodel --shift -75 --extsize 150 -B --keep-dup all --SPMR --outdir "${4}"/macs3 2> "${4}"/macs3/"${n:0:-4}".log" >> "${4}"/macs3_cmds
done
parallel -j $cores < "${4}"/macs3_cmds

# You want to go from .bedgraph to .bw and you can do this through kent utilities bedgraphToBigWig!
export PATH=/programs/kentUtils/bin:$PATH

# Now use a for loop to iterate through all the .bdg files and turn them into bigwigs!
> "${4}"/bw_cmds
for f in "${4}"/macs3/*_treat_pileup.bdg; do echo "bedGraphToBigWig $f /workdir/jdc397/1_currentWork/8_teCactus/0_hscData/mm10_sizes.txt "${f:0:-17}".bw" >> "${4}"/bw_cmds; done
parallel -j $cores < "${4}"/bw_cmds

# Finally, remove all the "control_lambda.bdg" and ".xls" files since you don't need nor want them! Compress the other .bdg!
rm "${4}"/macs3/*_control_lambda.bdg "${4}"/macs3/*.xls
> "${4}"/pigz_cmds
for f in $(ls "${4}"/macs3/*.bdg); do echo "pigz -p "${cores2}" $f" >> "${4}"/pigz_cmds; done
parallel -j $cores < "${4}"/pigz_cmds


######################################
### Shuffle for TE enrichment! 
######################################

# With the peaks generated, shuffle the summit calls for each peak using the power of your TEshuffle script, which should be broadly applicable
# to all datasets. Woot woot
module load R/4.1.3-r9
Rscript --vanilla TEshuffle.r "${4}"/macs3 "${4}"/shuffleResults *_peaks.narrowPeak 48

