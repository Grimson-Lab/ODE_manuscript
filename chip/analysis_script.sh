# !/bin/bash

# This script encapsulates all of your ChIP-seq pipelines for single and paired-end data using a manually added SRA Run Selector list of 
# data/labels. Use the first and second parameters to specify the specific "SRAfilt" R script, the input to that script, and the output 
# directory for everything!
# 
# Expected script call (from "ODE_manuscript directory")
# sh analysis_script [output directory]

### NOTE: You went through numerous papers and datasets to select the ones you ended up using. Regardless, their source can be found near the 
### end of the data list in the file "supData_1.txt" in the main directory

######################################
### Download the Data
######################################

# State the current task
echo "#######################"
echo "Downloading raw data..."
echo "#######################"

# Since you manually compiled the data file, generate the commands to get all of the files downloaded using awk and fasterq-dump
> "${1}"/dump_cmds

# Start with all the single-end read data
awk -v a=$1 '{FS="\t"; OFS=""; if ( NR != 1 && $3 == "Single" ) { print "fasterq-dump ", $1 " -o " a "/" $6 ".fq; pigz -p 4 " a "/" $6 ".fq" } }' "${1}"/SRA.tab > "${1}"/dump_cmds

# Now append all the paired-end read data
awk -v a=$1 '{FS="\t"; OFS=""; if ( NR != 1 && $3 == "Paired" ) { print "fasterq-dump ", $1 " -S; pigz -c -p 4 " $1 "_1.fastq > " a "/" $6 "a.fq.gz; pigz -c -p 4 " $1 "_2.fastq > " a "/" $6 "b.fq.gz; rm " $1 "_1.fastq " $1 "_2.fastq" } }' < "${1}"/SRA.tab >> "${1}"/dump_cmds

# Run the above commands to download the data!
parallel -j 6 < "${1}"/dump_cmds

# State the current task
echo "################################"
echo "Finished downloading raw data..."
echo "################################"



######################################
### Trim the Data (with Parallel!)
######################################

# Make the output directories
# mkdir "${1}"/trim; mkdir "${1}"/trim/fastqc

# Make sure to add the most recent cutadapt path which uses python3.0+ and allows for multi-threading!
export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages
export PATH=/programs/cutadapt-4.1/bin:$PATH

# State the current task
echo "########################"
echo "Starting the trimming..."
echo "########################"

# Write the trim commands to a file and then run those commands using parallel. Check to make sure the length is approximately 25% of the total
# read size (usually given as full read size, not the size of each mate), using the power of variables and math starting with single-end data
awk -v a=$1 '{FS="\t"; OFS=""; if ( NR != 1 && $3 == "Single" ) { print "trim_galore -j 4 --stringency 4 --length " int($2 / 4 + 0.5) " -o " a "/trim " a "/" $6 ".fq.gz --fastqc_args \"--outdir " a "/trim/fastqc/\" 2>/dev/null"} }' < "${1}"/SRA.tab > "${1}"/trim_cmds

# Now for the paired-end data
awk -v a=$1 '{FS="\t"; OFS=""; if ( NR != 1 && $3 == "Paired" ) { print "trim_galore -j 4 --stringency 4 --length " int($2 / 8 + 0.5) " --paired -o " a "/trim " a "/" $6 "a.fq.gz " a "/" $6 "b.fq.gz --fastqc_args \"--outdir " a "/trim/fastqc/\" 2>/dev/null"} }' < "${1}"/SRA.tab >> "${1}"/trim_cmds

# Run the above commands!
parallel -j 6 < "${1}"/trim_cmds

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
# mkdir "${1}"/STAR

# State the current task
echo "#####################"
echo "Starting alignment..."
echo "#####################"

# Populate the command file for single-end and then paired-end data and then align using parallel 
> "${1}"/star_cmds
for f in $(ls "${1}"/trim/*trimmed*.gz | sed 's/_trimmed.fq.gz//' | uniq)
do 
	n=`echo $f | sed 's|.*/||g'`
	echo "STAR --genomeDir STAR_noGene_mm10 --outFileNamePrefix "${1}"/STAR/"${n}" --alignMatesGapMax 2000 --outSAMtype BAM SortedByCoordinate --outFilterMatchNminOverLread .4 --seedSearchStartLmax 15 --readFilesCommand zcat --outFilterMultimapNmax 500 --outSAMmultNmax 1 --outMultimapperOrder Random --alignIntronMax 1 --readFilesIn "${f}"_trimmed.fq.gz --runMode alignReads --alignEndsType Local --runThreadN 4" >> "${1}"/star_cmds
done

for f in $(ls "${1}"/trim/*val*.gz | sed 's/[a-b]_val_[1-2].fq.gz//' | uniq)
do 
	n=`echo $f | sed 's|.*/||g'`
	echo "STAR --genomeDir STAR_noGene_mm10 --outFileNamePrefix "${1}"/STAR/"${n}" --alignMatesGapMax 2000 --outSAMtype BAM SortedByCoordinate --outFilterMatchNminOverLread .4 --seedSearchStartLmax 15 --readFilesCommand zcat --outFilterMultimapNmax 500 --outSAMmultNmax 1 --outMultimapperOrder Random --alignIntronMax 1 --readFilesIn "${f}"a_val_1.fq.gz "${f}"b_val_2.fq.gz --runMode alignReads --alignEndsType Local --runThreadN 4" >> "${1}"/star_cmds
done
parallel -j 6 < "${1}"/star_cmds

# State the current task
echo "#####################"
echo "Finished alignment..."
echo "#####################"


# State the current task
echo "##############################"
echo "Starting merging & deduping..."
echo "##############################"

# Create a file of all the single-end and paired-end data for below, since they aren't easily distinguishable with file names anymore
sed '1d' "${1}"/SRA.tab | awk '{OFS="\t"; if ($3 == "Single") print $6}' | sed 's/_[1-9]$//g' | sort | uniq > "${1}"/single_data
sed '1d' "${1}"/SRA.tab | awk '{OFS="\t"; if ($3 == "Paired") print $6}' | sed 's/_[1-9]$//g' | sort | uniq > "${1}"/paired_data

# Now that you have the aligned data, merge all of the replicates (if present) and then run Picard to remove the PCR/Optical duplicates as 
# well as filter out all of the reads that align to the gross blacklist regions and chrM/Y
> "${1}"/picard_cmds
while read f
do
	g=`echo ""${1}"/STAR/"${f}""`
	if [[ `ls "${g}"_[1-9]*Aligned* | wc -l` -gt 1 ]]
	then 
		echo "samtools merge -@4 -o "$g"_merge.bam $(ls "${g}"_[1-9]*.out.bam | tr '\n' ' '); samtools view -@4 -h -q 10 "$g"_merge.bam | awk '{if (\$3 != \"chrM\" && \$3 != \"chrY\") print}' | samtools view -b -@4 | bedtools intersect -v -wa -a stdin -b comb_badRegions.bed > "$g"_tmp.bam; java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates I="$g"_tmp.bam ASSUME_SORT_ORDER=coordinate O="$g".dedup.bam M="$g".metrics REMOVE_DUPLICATES=true" >> "${1}"/picard_cmds
	else 
		echo "samtools view -@4 -h -q 10 $(ls "${g}"_[1-9]*.out.bam) | awk '{if (\$3 != \"chrM\" && \$3 != \"chrY\") print}' | samtools view -b -@4 | bedtools intersect -v -wa -a stdin -b comb_badRegions.bed > "$g"_tmp.bam; java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates I="$g"_tmp.bam ASSUME_SORT_ORDER=coordinate O="$g".dedup.bam M="$g".metrics REMOVE_DUPLICATES=true" >> "${1}"/picard_cmds
	fi

done < "${1}"/single_data

# Since the only paired data you have has replicates, you don't need the if/else as above
while read f
do
	g=`echo ""${1}"/STAR/"${f}""`
	echo "samtools merge -@4 -o "$g"_merge.bam $(ls "${g}"*.out.bam | tr '\n' ' '); samtools view -@4 -hf 2 -q 10 "$g"_merge.bam | awk '{if (\$3 != \"chrM\" && \$3 != \"chrY\") print}' | samtools view -b -@4 | bedtools intersect -v -wa -a stdin -b comb_badRegions.bed > "$g"_tmp.bam; java -jar /programs/picard-tools-2.19.2/picard.jar MarkDuplicates I="$g"_tmp.bam ASSUME_SORT_ORDER=coordinate O="$g".dedup.bam M="$g".metrics REMOVE_DUPLICATES=true" >> "${1}"/picard_cmds

done < "${1}"/paired_data

parallel -j 6 < "${1}"/picard_cmds

# State the current task
echo "##############################"
echo "Finished merging & deduping..."
echo "##############################"

# Remove the temporary files 
rm "${1}"/STAR/*_tmp.bam "${1}"/STAR/*Aligned.sortedByCoord.out.bam "${1}"/STAR/*_merge.bam

# Remove all the original .fq.gz files and the intermediate .bam files you no longer need
rm "${1}"/*.fq.gz; rm "${1}"/STAR/*.tab



######################################
### Macs3 the Data 
######################################

# Use the power of MACS3 to call peaks using the relevant input files (if they exist)
source ~/Macs3Env/bin/activate 

# State the current task
echo "################################"
echo "Starting peak calling w/MACS3..."
echo "################################"

# Make the output directory
# mkdir macs3 

# Fill the commands file with commands.  
> "${1}"/macs3_cmds
for f in $(ls "${1}"/STAR/*.bam | grep -v "input" | sed 's/.dedup.bam//g' | uniq)
do 
	if [[ -f "${f}"_input.dedup.bam ]]
	then 
		echo "macs3 callpeak -t "${f}".dedup.bam -c "${f}"_input.dedup.bam -f AUTO -g mm -n "${f:10}" -q 0.05 -B --keep-dup all --call-summits --SPMR --outdir "${1}"/macs3 2> "${1}"/macs3/"${f:10}".log" >> "${1}"/macs3_cmds
	else
		echo "macs3 callpeak -t "${f}".dedup.bam -f AUTO -g mm -n "${f:10}" -q 0.05 -B --keep-dup all --call-summits --SPMR --outdir "${1}"/macs3 2> "${1}"/macs3/"${f:10}".log" >> "${1}"/macs3_cmds
	fi
	
done
parallel -j 6 < "${1}"/macs3_cmds

# State the current task
echo "################################"
echo "Finished peak calling w/MACS3..."
echo "################################"


# State the current task
echo "############################"
echo "Starting file compression..."
echo "############################"

# You want to go from .bedgraph to .bw and you can do this through kent utilities bedgraphToBigWig!
export PATH=/programs/kentUtils/bin:$PATH

# Iterate through all the sample .bdg files and turn them into bigwigs!
> "${1}"/bw_cmds
for f in "${1}"/macs3/*_treat_pileup.bdg; do echo "bedGraphToBigWig $f mm10_sizes.txt "${f:0:-17}".bw" >> "${1}"/bw_cmds; done
parallel -j 6 < "${1}"/bw_cmds

# Compress all of the .narrowPeak files
gzip "${1}"/macs3/*.narrowPeak

# Finally, remove all the superfluous files!
rm "${1}"/macs3/*.bdg "${1}"/macs3/*.xls "${1}"/macs3/*model.r

# State the current task
echo "############################"
echo "Finished file compression..."
echo "############################"
