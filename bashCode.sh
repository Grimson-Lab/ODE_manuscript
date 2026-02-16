### This code will download the mouse .fa file and generate a STAR (v2.7.5a) genome without gene annotations for use in aligning 
### RNA-seq, ATAC-seq, and ChIP-seq data

# First, download the mm10 genome in fasta sequnce and then unzip it
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# You don't want any of the unknown or random contigs, so remove them using samtools 
samtools faidx mm10.fa
grep ">" mm10.fa | egrep -v 'chrU|_random' | sed 's/^>//g' | xargs samtools faidx mm10.fa > clean_mm10.fa

# Rename the filtered file to overwrite the original
rm mm10.fa.fai; mv clean_mm10.fa mm10.fa

### NOTE: On our server, we use the below code to set the right version of STAR as default for the session
# export PATH=/programs/STAR-2.7.5a/bin/Linux_x86_64_static:$PATH 

# Now, use STAR to generate the genome index file for all alignment needs throughout the analysis
STAR --runMode genomeGenerate --runThreadN 24 --genomeDir STAR_noGene_mm10 --genomeFastaFiles mm10.fa

# Build a bowtie2 index for identifying regions similar to chrM across the genome
bowtie2-build -f --threads 24 mm10.fa mm10; mkdir bowtie2_mm10; mv *.bt2 bowtie2_mm10/.



## For the alignments, you want to filter out reads that fall within regions of high input signal, mitochondrial homology, or satellites. 
## As for the Mitochondrial homology, generate artificial reads from the mitochondrial chromosome to see what sequence is similar

# Download just the mouse mitochondrial genome to generate reads from just that 'chromosome'
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chrM.fa.gz; gunzip chrM.fa.gz

# Now generate an artificial .fastq file so you can align them back to the genome! 
java -jar ~/ArtificialFastqGenerator/ArtificialFastqGenerator.jar -R chrM.fa -O chrM -S ">chrM" 

# Use Bowtie2 for more flexible mapping (vs STAR) of the above reads  
bowtie2 --very-sensitive-local -a --no-unal -x bowtie2_mm10/mm10 -1 chrM.1.fastq -2 chrM.2.fastq -p 8 2>chrM.log | samtools -@8 view -bS > chrM.bam

# Now use bamtobed to identify all the regions in the genome where homology with the mitochondria exists
bedtools bamtobed -i chrM.bam | grep -v "chrM" | sort -k1,1 -k2,2n | bedtools merge -d 200 > chrM_badRegions.bed

# Download the ENCODE blacklist for mm10 from this paper: 
# https://www.nature.com/articles/s41598-019-45839-z
wget https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz; gunzip ENCFF547MET.bed.gz

# Now combine the mitochondrial bad regions with the above blacklist and annotated satellite regions
cat chrM_badRegions.bed ENCFF547MET.bed <(zcat rm/mm10_repeatMasker.txt | grep "Satellite" | cut -f1-3) | sort -k1,1 -k2,2n | bedtools merge -d 100 > comb_badRegions.bed


# Download the list of gene annotations for mm10 to use throughout the analysis for RNA-seq and ABC model: RefSeq Curated
# 
# You clicked the following buttons: Clade: Mammal, Genome: Mouse, Assembly: mm10, Group: Genes and Gene Predicitions, Track: NCBI RefSeq, 
# Table: RefSeq Curated, Output format: GTF and then you downloaded the file (compressed) as "mm10.refseqCurated.gtf.gz"
# 
# gunzip mm10.refseqCurated.gtf.gz

# Download the .txt file containing all of the name2 information for the above gene annotations from UCSC 
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqCurated.txt.gz 


### Create a chromosome file without chrM and chrY
egrep -v 'chrM|chrY' STAR_noGene_mm10/chrNameLength.txt > mm10_sizes.txt



#############################
### Run Analysis Scripts! ###
#############################

### Run RepeatMasker on mm10 (takes a few days...)

sh rm/analysis_script.sh rm


### Download and analyze the ImmGen ATAC-seq

sh atac/immgen/analysis_script.sh SRAfilt.r SraRunTable_immgenATAC.tab SRA.tab atac/immgen


### Download and analyze the ENCODE ATAC-seq 

sh atac/encode/analysis_script.sh SRAfilt.r SraRunTable_encodeATAC.csv SRA.tab atac/encode


### Download and analyze the ENCODE DNase-seq

sh dnase/analysis_script.sh


### Download and analyze the ImmGen RNA-seq

sh rna/analysis_script.sh SRAfilt.r SraRunTable_immgenRNA.tab SRA.tab rna


### Download and analyze the transcription factor ChIP-seq

sh chip/analysis_script.sh chip



### NOTE: With the above data downloaded and analyzed, run all of the necessary bash code to get the files needed for each figure's panels

### NOTE: The below code is written with the assumption that all files generated from the above bash scripts are present or otherwise provided
### on the github page. 


################
### Figure 1 ###
################

# Generate a file containing each of the ImmGen datasets to iterate over for downstream analyses. Yeah! 
sed '1d' atac/immgen/TEshuffle_results.txt | cut -f13 | sed 's/_onlyEnhancer.txt//g' | sort | uniq > atac/immgen.iter

# Let's get that unified peak set! You used the median peak size to generate sample-specific peak sizes off of the summit using code from
# here: https://stackoverflow.com/questions/6166375/median-of-column-with-awk and you also removed all sample peaks with size >2kb at first, 
# and then removed all unified peaks >10kb to minimize dealing with super long regions that would mess up downstream analyses 
> tmp
for f in $(ls atac/immgen/macs3/*.narrowPeak); do val=`awk '{if ($3 - $2 <= 2000) print $3 - $2}' $f | sort -n | awk '{ a[i++] = $1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print int(a[x-1]/2) }'`; awk '{OFS="\t"; if ($3 - $2 <= 2000) print $1, $2 + $10, $2 + $10 + 1}' $f | bedtools slop -b $val -i stdin -g mm10_sizes.txt >> tmp; done
sort -k1,1 -k2,2n -S 10% --parallel=24 tmp | bedtools merge | awk '{OFS="\t"; if ($3 - $2 <= 10000) print}' | awk '{OFS="\t"; print $0, "peak_" ++count}' > atac/summit_unif_peaks.bed
rm tmp

# Generate the above file with rownames
awk '{OFS="\t"; print $1 ":" $2 "-" $3 "|" $4, $0}' atac/summit_unif_peaks.bed > atac/full_summit_unif_peaks.bed

# Automate the amount of fragments that map to each peak for each sample & replicate into a single file! Yay!
> atac/int_cmds
for f in $(ls atac/immgen/STAR/*.bed)
do
	n=`echo $f | sed 's|.*/||g' | sed 's|.bed||g'`
	echo "sort -k1,1 -k2,2n $f | intersectBed -sorted -wa -c -a atac/summit_unif_peaks.bed -b stdin | cut -f5 > $n" >> atac/int_cmds 
done 
parallel -j 30 < atac/int_cmds

# Append all the intersected files together. Yeah for speed
cat atac/summit_unif_peaks.bed > tmp
while read f; do paste tmp $f > tmp2; cat tmp2 > tmp; rm $f; done < atac/immgen.iter
cat tmp > atac/summit_unif_peaks_counts.txt; rm tmp tmp2

# With the peak files created, it's time to identify which TE(s) are associated with which peaks based on the summit overlap across each 
# dataset so that you can assign them later on with downstream analyses
> atac/teint_cmds
while read f; do echo "awk '{OFS=\"\\t\"; print \$1, \$2 + \$10, \$2 + \$10 + 1}' atac/immgen/macs3/"${f}"_peaks.narrowPeak | sort -k1,1 -k2,2n | intersectBed -sorted -wa -wb -a stdin -b rm/mm10_repeatMasker.txt.gz | awk '{OFS=\"\\t\"; print \$0, \""$f"\"}' > "$f"" >> atac/teint_cmds; done < atac/immgen.iter
parallel -j 32 < atac/teint_cmds

# Read through the above files, combine them together into a single file, and then remove all the read-in files for space
> atac/te.bound
while read f; do cat $f >> atac/te.bound; rm $f; done < atac/immgen.iter

# Create a list of all of the TEs that have a summit overlap them across any of the above datasets. Yeah! 
awk '{OFS="\t"; print $4, $5, $6, $4 ":" $5 "-" $6 "|" $7, $8, $9}' atac/te.bound | sort -k1,1 -k2,2n -S 10% --parallel=8 | uniq > atac/accessible.tes.atac.txt

# Finally, generate the list of TE-associated peaks based on the above list of accessible elements. A peak is TE-derived if any of the summits 
# from a dataset overlap a TE that makes up the peak. Yeah
cut -f1-11 atac/te.bound | sort -k1,1 -k2,2n -S 10% --parallel=16 | uniq | intersectBed -sorted -wa -wb -a atac/summit_unif_peaks.bed -b stdin | cut -f1-4,8-15 | uniq > atac/peak.te.int



################
### Figure 2 ###
################

## B

# Add HOMER to your path (can be downloaded and initialized from http://homer.ucsd.edu/homer/configureHomer.pl)
export PATH=/workdir/jdc397/homer/bin:$PATH

# Iterate through the maximal number of groups based on whether there are elements within them. Yeah
for i in {1..8}
do
	# Start with ORR1E
	awk -v a=$i '{FS="\t"; OFS="\t"; if ($7 == a && $4 ~ "1E") print}' atac/rObject_ode.txt | sort -k1,1 -k2,2n --parallel=6 -S 10% | uniq | intersectBed -sorted -wa -f .95 -a rm/mm10_repeatMasker.txt.gz -b stdin | awk '{OFS="\t"; print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' | fastaFromBed -s -nameOnly -fi mm10.fa -bed stdin | sed 's/(.*//' > tmp
	for j in {1..8}
	do
		awk -v a=$j '{FS="\t"; OFS="\t"; if ($7 == a && $4 ~ "1E") print}' atac/rObject_ode.txt | sort -k1,1 -k2,2n --parallel=6 -S 10% | uniq | intersectBed -sorted -wa -f .95 -a rm/mm10_repeatMasker.txt.gz -b stdin | awk '{OFS="\t"; print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' | fastaFromBed -s -nameOnly -fi mm10.fa -bed stdin | sed 's/(.*//' > tmp2
		if [[ `cat tmp | wc -l` != 0 && `cat tmp2 | wc -l` != 0 && $j -ne $i ]]
		then
			homer2 known -i tmp -b tmp2 -o homer/mouse/orr1e_km_$i.vs.$j.txt -m homer/homerFilt.motifs -p 24
		else
			echo "Skipping group $i & $j for ORR1E"
		fi
	done

	# Now for ORR1D2
	awk -v a=$i '{FS="\t"; OFS="\t"; if ($7 == a && $4 ~ "1D2") print}' atac/rObject_ode.txt | sort -k1,1 -k2,2n --parallel=6 -S 10% | uniq | intersectBed -sorted -wa -f .95 -a rm/mm10_repeatMasker.txt.gz -b stdin | awk '{OFS="\t"; print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' | fastaFromBed -s -nameOnly -fi mm10.fa -bed stdin | sed 's/(.*//' > tmp
	for j in {1..6}
	do
		awk -v a=$j '{FS="\t"; OFS="\t"; if ($7 == a && $4 ~ "1D2") print}' atac/rObject_ode.txt | sort -k1,1 -k2,2n --parallel=6 -S 10% | uniq | intersectBed -sorted -wa -f .95 -a rm/mm10_repeatMasker.txt.gz -b stdin | awk '{OFS="\t"; print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' | fastaFromBed -s -nameOnly -fi mm10.fa -bed stdin | sed 's/(.*//' > tmp2
		if [[ `cat tmp | wc -l` != 0 && `cat tmp2 | wc -l` != 0 && $j -ne $i ]]
		then
			homer2 known -i tmp -b tmp2 -o homer/mouse/orr1d2_km_$i.vs.$j.txt -m homer/homerFilt.motifs -p 24
		else
			echo "Skipping group $i & $j for ORR1D2"
		fi
	done
done
rm tmp tmp2


# Run TOMTOM from the MEME-suite on the vertebrate HOMER motif file (turned into MEME format using the universalmotif package in R)
/programs/meme-5.5.2/bin/tomtom -dist kullback -motif-pseudo 0.1 -text -min-overlap 1 homer/homerFilt.meme homer/homerFilt.meme > homer/tomtom.homer.all.txt

# Run the python script to take the above and cluster each of the motifs based on positional-weight matrix
python homer/motifCluster.py


## C

# Use the file with all accessible ODEs to create the foundation for the below heatmaps
for f in {1..8}
do 
	if [[ f -eq 1 ]]
	then
		awk -v a=$f '{OFS="\t"; if ($4 == "ORR1E" && $7 == a) print}' atac/rObject_ode.txt | cut -f1-3 | uniq > tmp
	else
		awk -v a=$f '{OFS="\t"; if ($4 == "ORR1E" && $7 == a) print}' atac/rObject_ode.txt | cut -f1-3 | uniq >> tmp
	fi 
	echo "#${f}" >> tmp
done
cat tmp > chip/orr1e_heatmap_base.bed
rm tmp 

# Prep the paths for deepTools (which is what was done on our server. You may/may not need to do your own version of below)!
export PATH=/programs/deeptools-3.5.5/bin:$PATH
export PYTHONPATH=/programs/deeptools-3.5.5/lib64/python3.9/site-packages:/programs/deeptools-3.5.5/lib/python3.9/site-packages
export TMPDIR=/workdir/$USER/tmp

# Generate the heatmaps with default parameters to guide the scaling between different heatmaps
computeMatrix scale-regions -S chip/macs3/prog_pu1.bw -R chip/orr1e_heatmap_base.bed -o orr1e_preHeatmap_pc --outFileNameMatrix orr1e_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1e_preHeatmap_pc --outFileSortedRegions prog_pu1.txt -o "orr1e_prog_pu1.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/matureB_pax5.bw -R chip/orr1e_heatmap_base.bed -o orr1e_preHeatmap_pc --outFileNameMatrix orr1e_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1e_preHeatmap_pc --outFileSortedRegions matureB_pax5.txt -o "orr1e_matureB_pax5.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/cd8T_runx1.bw -R chip/orr1e_heatmap_base.bed -o orr1e_preHeatmap_pc --outFileNameMatrix orr1e_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1e_preHeatmap_pc --outFileSortedRegions cd8T_runx1.txt -o "orr1e_cd8T_runx1.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/th17_rorg.bw -R chip/orr1e_heatmap_base.bed -o orr1e_preHeatmap_pc --outFileNameMatrix orr1e_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1e_preHeatmap_pc --outFileSortedRegions th17_rorg.txt -o "orr1e_th17_rorg.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/bmdm_pu1.bw -R chip/orr1e_heatmap_base.bed -o orr1e_preHeatmap_pc --outFileNameMatrix orr1e_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1e_preHeatmap_pc --outFileSortedRegions bmdm_pu1.txt -o "orr1e_bmdm_pu1.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/cd8DC_irf8.bw -R chip/orr1e_heatmap_base.bed -o orr1e_preHeatmap_pc --outFileNameMatrix orr1e_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1e_preHeatmap_pc --outFileSortedRegions cd8DC_irf8.txt -o "orr1e_cd8DC_irf8.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

# Combine all the pdfs you made above into a single pdf for simplicity using the power of pdfunite 
# pdfunite $(ls orr1e*.pdf) comb_tf_heatmap.pdf

# Remove the above .pdfs since you only need the sorted regions
rm $(ls orr1e*.pdf)

# Combine the sorted rows based on the respective ChIP-seq dataset
awk '{OFS="\t"; if ($13 == 1 || $1 == "#chrom") print}' prog_pu1.txt > comb_chip_orr1e_sort.txt
awk '{OFS="\t"; if ($13 == 2) print}' matureB_pax5.txt >> comb_chip_orr1e_sort.txt
awk '{OFS="\t"; if ($13 == 3 || $13 == 4) print}' cd8T_runx1.txt >> comb_chip_orr1e_sort.txt
awk '{OFS="\t"; if ($13 == 5) print}' th17_rorg.txt >> comb_chip_orr1e_sort.txt
awk '{OFS="\t"; if ($13 == 6) print}' bmdm_pu1.txt >> comb_chip_orr1e_sort.txt
awk '{OFS="\t"; if ($13 == 7 || $13 == 8) print}' cd8DC_irf8.txt >> comb_chip_orr1e_sort.txt

computeMatrix scale-regions -R comb_chip_orr1e_sort.txt -S chip/macs3/prog_pu1.bw chip/macs3/matureB_pax5.bw chip/macs3/cd8T_runx1.bw chip/macs3/th17_rorg.bw chip/macs3/bmdm_pu1.bw chip/macs3/cd8DC_irf8.bw -b 1000 -a 1000 -o orr1e_preHeatmap_pc --outFileNameMatrix orr1e_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1e_preHeatmap_pc --sortRegions no --zMin 0 --zMax 1.5 1 1.5 0.3 3 1 --outFileSortedRegions chip/comb_chip_orr1e.txt -o "figures/fig2c_orr1eIndividualScale_heatmap.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

# Now remove the temporary files
rm orr1e_preHeatmap_pc orr1e_preHeatmap_pc.mat prog_pu1.txt matureB_pax5.txt cd8T_runx1.txt th17_rorg.txt bmdm_pu1.txt cd8DC_irf8.txt



################
### Figure 3 ###
################

## A

# With the above complete, it's time to generate an upset plot for figure 3! To do this, you need to get a table with all of the overlaps 
# between ODEs and the datasets that you chose for the figure. You manually input those into a file named "samps_fig3" 
while read f
do 
	name=`echo $f | sed 's|.*/||g' | sed 's/_peaks.*//g'`
	sort -k1,1 -k2,2n atac/rObject_ode.txt | sed '1d' | cut -f1-4,7 | uniq | intersectBed -sorted -wa -c -a stdin -b $f | cut -f6 > $name
done < chip/samps_fig3

sed '1d' atac/rObject_ode.txt | sort -k1,1 -k2,2n | cut -f1-5,7 | uniq > tmp; > chip/figure3_chip_names
while read f
do 
	name=`echo $f | sed 's|.*/||g' | sed 's/_peaks.*//g'`; echo $name >> chip/figure3_chip_names
	paste tmp $name > tmp2; cat tmp2 > tmp; rm $name
done < chip/samps_fig3
cat tmp > chip/figure3_ode_int_chip_vals.txt; rm tmp tmp2


## F

# Generate a fasta file for the cell type-specific ODE elements 
sort -k1,1 -k2,2n --parallel=6 -S 10% atac/rObject_ode.txt | uniq | intersectBed -sorted -wa -f .95 -a rm/mm10_repeatMasker.txt.gz -b stdin | awk '{OFS="\t"; print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' | fastaFromBed -s -nameOnly -fi mm10.fa -bed stdin | sed 's/(.*//' > rm/ode.fa

# Identify the ODE elements which are ~full length and could be used as a background 
awk '{OFS="\t"; if ($4 == "ORR1E" && $3 - $2 >= 341 && $3 - $2 <= 377) print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' rm/mm10_repeatMasker.txt.gz | intersectBed -sorted -v -wa -a stdin -b atac/accessible.tes.atac.txt > rm/fl_ode.txt
awk '{OFS="\t"; if ($4 == "ORR1D2" && $3 - $2 >= 351 && $3 - $2 <= 389) print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' rm/mm10_repeatMasker.txt.gz | intersectBed -sorted -v -wa -a stdin -b atac/accessible.tes.atac.txt >> rm/fl_ode.txt

# Turn the above full-length non-accessible elements into a fasta file
fastaFromBed -s -nameOnly -fi mm10.fa -bed rm/fl_ode.txt | sed 's/(.*//' > rm/fl_ode.fa


# Download phastCons scores from the Glire subset of the 60way comparison, which includes 8 total species. You believe ORR1 should
# be present in all of these species, except possibly Pika (Ochotona princeps) since it's split the farthest out from the rest of the species.
rsync rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons60wayGlire.bw rm/.

# Turn the above into a .bedgraph for intersection purposes
export PATH=/programs/kentUtils/bin:$PATH
bigWigToBedGraph rm/mm10.60way.phastCons60wayGlire.bw rm/mm10.60way.phastCons60wayGlire.bedgraph 

# Generate the intersected .bedgraph files for glire phastCons scores across the ORR1 elements for loading into R 
awk '{OFS="\t"; if ($4 == "ORR1E" || $4 == "ORR1D2") print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' rm/mm10_repeatMasker.txt.gz | intersectBed -sorted -wa -wb -a stdin -b rm/mm10.60way.phastCons60wayGlire.bedgraph | pigz -c -p 24 > orthology/ode.phastCons60way.int.txt.gz



################
### Figure 6 ###
################

## B

rsync rsync://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.fa.gz orthology/.; gunzip orthology/rn7.fa.gz
makeblastdb -in orthology/rn7.fa -blastdb_version 5 -title "rn7" -dbtype nucl &

# First, slop the output, this time by 2kb on either side. Then, blast it!
slopBed -b 1000 -i <( sed '1d' atac/rObject_ode.txt | sort -k1,1 -k2,2n --parallel=6 -S 10% | intersectBed -sorted -wa -f .95 -a rm/mm10_repeatMasker.txt.gz -b stdin | uniq | awk '{OFS="\t"; print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4, $5, $6}' ) -g mm10_sizes.txt | sort -k1,1 -k2,2n > ode_s1000.bed 
fastaFromBed -nameOnly -s -fi mm10.fa -bed ode_s1000.bed | sed 's/(.*//' > ode_s1000.fa
blastn -num_threads 64 -evalue 1e-30 -qcov_hsp_perc 15 -outfmt "6 qseqid sseqid pident nident length mismatch qstart qend sstart send qlen evalue bitscore" -query ode_s1000.fa -db orthology/rn7.fa -out orthology/rn7_1000_blast_ode.txt

# Now that all of the blasts were run, compile them together into a single object for easy(ish) loading into R
awk '{OFS="\t"; if ( $10 > $9 ) { print $0, "+" } else { print $1, $2, $3, $4, $5, $6, $7, $8, $10, $9, $11, $12, $13, "-" } }' orthology/rn7_1000_blast_ode.txt | sort -k2,2 -k9,9n --parallel=32 -S 20% > tmp_sort
Rscript --vanilla scripts/blastMerge.r tmp_sort orthology/rn7_1000_blast_ode_output.txt

# Remove the tmp file and all other superfluous files
rm tmp_sort ode_s1000.bed ode_s1000.fa

# To prep for the alignments, get those going in R
awk '{OFS="\t"; print $2, $9, $10, $1, $3, $14, $7, $8}' orthology/rn7_1000_blast_ode_output.txt | sort -k1,1 -k2,2n > orthology/ode_blast_s1000.bed
fastaFromBed -nameOnly -s -fi orthology/rn7.fa -bed orthology/ode_blast_s1000.bed | sed 's/(.*//' > orthology/ode_blast_s1000.fa

### NOTE: R code found in Figure 6B generates the relevant .fa files that the below code uses and are provided by yours truly 

# Add HOMER to your path
export PATH=/workdir/jdc397/homer/bin:$PATH

# Iterate through the maximal number of groups based on whether there are elements within them. Yeah
for i in {1..8}
do
	# Start with ORR1E
	for j in {1..8}
	do
		if [[ $j -ne $i ]]
		then
			homer2 known -i orthology/ORR1E_blast_cpm.c_$i.fa -b orthology/ORR1E_blast_cpm.c_$j.fa -o homer/rat/orr1e_km_$i.vs.$j.txt -m homer/homerFilt.motifs -p 24
		else
			echo "Skipping group $i & $j for ORR1E"
		fi
	done
done


# Now for ORR1D2
for i in {1..6}
do
	for j in {1..6}
	do
		if [[ $j -ne $i ]]
		then
			homer2 known -i orthology/ORR1D2_blast_cpm.c_$i.fa -b orthology/ORR1D2_blast_cpm.c_$j.fa -o homer/rat/orr1d2_km_$i.vs.$j.txt -m homer/homerFilt.motifs -p 24
		else
			echo "Skipping group $i & $j for ORR1D2"
		fi
	done
done



#############################
### Supplemental Figure 2 ###
#############################

## E

# Get the list of ORR1D2s for the below heatmap
for f in {1..6}
do 
	if [[ f -eq 1 ]]
	then
		awk -v a=$f '{OFS="\t"; if ($4 == "ORR1D2" && $7 == a) print}' atac/rObject_ode.txt | cut -f1-3 | uniq > tmp
	else
		awk -v a=$f '{OFS="\t"; if ($4 == "ORR1D2" && $7 == a) print}' atac/rObject_ode.txt | cut -f1-3 | uniq >> tmp
	fi 
	echo "#${f}" >> tmp
done
cat tmp > chip/orr1d2_heatmap_base.bed
rm tmp

# Prep the paths for deepTools!
export PATH=/programs/deeptools-3.5.5/bin:$PATH
export PYTHONPATH=/programs/deeptools-3.5.5/lib64/python3.9/site-packages:/programs/deeptools-3.5.5/lib/python3.9/site-packages
export TMPDIR=/workdir/$USER/tmp

# Generate the heatmaps with default parameters to guide the scaling between different heatmaps
computeMatrix scale-regions -S chip/macs3/prog_pu1.bw -R chip/orr1d2_heatmap_base.bed -o orr1d2_preHeatmap_pc --outFileNameMatrix orr1d2_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1d2_preHeatmap_pc --outFileSortedRegions prog_pu1.txt -o "orr1d2_prog_pu1.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/matureB_pax5.bw -R chip/orr1d2_heatmap_base.bed -o orr1d2_preHeatmap_pc --outFileNameMatrix orr1d2_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1d2_preHeatmap_pc --outFileSortedRegions matureB_pax5.txt -o "orr1d2_matureB_pax5.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/cd8T_runx1.bw -R chip/orr1d2_heatmap_base.bed -o orr1d2_preHeatmap_pc --outFileNameMatrix orr1d2_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1d2_preHeatmap_pc --outFileSortedRegions cd8T_runx1.txt -o "orr1d2_cd8T_runx1.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/th17_rorg.bw -R chip/orr1d2_heatmap_base.bed -o orr1d2_preHeatmap_pc --outFileNameMatrix orr1d2_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1d2_preHeatmap_pc --outFileSortedRegions th17_rorg.txt -o "orr1d2_th17_rorg.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/bmdm_pu1.bw -R chip/orr1d2_heatmap_base.bed -o orr1d2_preHeatmap_pc --outFileNameMatrix orr1d2_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1d2_preHeatmap_pc --outFileSortedRegions bmdm_pu1.txt -o "orr1d2_bmdm_pu1.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

computeMatrix scale-regions -S chip/macs3/cd8DC_irf8.bw -R chip/orr1d2_heatmap_base.bed -o orr1d2_preHeatmap_pc --outFileNameMatrix orr1d2_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1d2_preHeatmap_pc --outFileSortedRegions cd8DC_irf8.txt -o "orr1d2_cd8DC_irf8.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

# Combine all the pdfs you made above into a single pdf for simplicity using the power of pdfunite 
# pdfunite $(ls orr1d2*.pdf) comb_tf_heatmap.pdf

# Remove the above .pdfs since you only need the sorted regions
rm $(ls orr1d2*.pdf)

# Combine the sorted rows based on the respective ChIP-seq dataset
awk '{OFS="\t"; if ($13 == 1 || $1 == "#chrom") print}' prog_pu1.txt > comb_chip_orr1d2_sort.txt
awk '{OFS="\t"; if ($13 == 2) print}' matureB_pax5.txt >> comb_chip_orr1d2_sort.txt
awk '{OFS="\t"; if ($13 == 3 || $13 == 4) print}' cd8T_runx1.txt >> comb_chip_orr1d2_sort.txt
awk '{OFS="\t"; if ($13 == 5) print}' th17_rorg.txt >> comb_chip_orr1d2_sort.txt
awk '{OFS="\t"; if ($13 == 6) print}' bmdm_pu1.txt >> comb_chip_orr1d2_sort.txt
awk '{OFS="\t"; if ($13 == 7 || $13 == 8) print}' cd8DC_irf8.txt >> comb_chip_orr1d2_sort.txt

computeMatrix scale-regions -R comb_chip_orr1d2_sort.txt -S chip/macs3/prog_pu1.bw chip/macs3/matureB_pax5.bw chip/macs3/cd8T_runx1.bw chip/macs3/th17_rorg.bw chip/macs3/bmdm_pu1.bw chip/macs3/cd8DC_irf8.bw -b 1000 -a 1000 -o orr1d2_preHeatmap_pc --outFileNameMatrix orr1d2_preHeatmap_pc.mat -p 24
plotHeatmap -m orr1d2_preHeatmap_pc --sortRegions no --zMin 0 --zMax 1.5 1 1.5 0.3 3 1 --outFileSortedRegions chip/comb_chip_orr1d2.txt -o "figures/fig2c_orr1d2IndividualScale_heatmap.pdf" --colorMap "Purples" --startLabel "5'" --endLabel "3'" --legendLocation "none" --whatToShow "heatmap and colorbar"

# Now remove the temporary files
rm orr1d2_preHeatmap_pc orr1d2_preHeatmap_pc.mat prog_pu1.txt matureB_pax5.txt cd8T_runx1.txt th17_rorg.txt bmdm_pu1.txt cd8DC_irf8.txt











### Updating is still in progress...

## Determine enrichment for the overlap of ABC-called CRE-gene linkages with the ATAC-seq and RNA-seq correlation linkages

# zcat /workdir/jdc397/microC/abc/immgen_cd8_rss/Predictions/EnhancerPredictionsAllPutative.tsv.gz | awk '{OFS="\t"; if ($5 != "promoter" && $28 >= 0.025) print $1, $2, $3, $4, $5, $11, $12, $16, $28}' > rss_filt_abc_links.txt

# cut -f1-4 atac/summit_unif_peaks_10k.txt | intersectBed -wa -wb -a stdin -b <( sed '1d' rss_filt_abc_links.txt | cut -f1-3,6,7 ) | awk '{OFS="\t"; print $1, $9, $9, $8, $5, $6, $7, $1, $2, $3, $4}' | intersectBed -wa -wb -a stdin -b refSC_mm10_tss.bed | awk '{OFS="\t"; print $8, $9, $10, $11, int(($14 + $13 + 1) / 2), $15}' | sort -u | sort -k1,1 -k2,2n -k5,5 > rss_abc_links_int_base_shuffle.txt

module load R/4.1.3-r9
Rscript --vanilla linkageShuffle.R rss_abc_links_int_base_shuffle.txt 64

# Code you used to go from the ABC model predictions to a filtered list of interactions
zcat /workdir/jdc397/microC/abc/immgen_cd8_rss/Predictions/EnhancerPredictionsAllPutative.tsv.gz | awk '{OFS="\t"; if ($5 != "promoter" && $28 >= 0.025) print $1, $2, $3, $4, $5, $11, $12, $16, $28}' > rss_filt_abc_links.txt

# To load in all the filtered contacts regardless of ODE intersection, get the below file generated
zcat microc/rss_filt_abc_links.txt.gz | intersectBed -wa -wb -a atac/summit_unif_peaks.bed -b stdin | awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $10, $12, $13}' | pigz -p 6 -c > microc/rss_abc_peak_int_all.txt.gz







