### Friday December 26th, 2025
# This is a compendium of the R code used to perform any supplemental analyses, such as measuring enrichment for the CRE-gene linkages between
# the ABC model and the ATAC-seq and RNA-seq correlation method. 



## Measuring enrichment of ABC-linkages in the ATAC-RNA correlation linkages

# Add lbzip2 to your path and start R version 4.1.3
# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; R --vanilla

# Load in all the packages you need
suppressPackageStartupMessages({ 
	library(fastSave); library(regioneR); library(data.table); setDTthreads(8)
})

# First, take the filtered correlations and output all of the ATAC-RNA linkages as a .bed file
load.lbzip2(file = "atac/rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Output a file that contains each of the linkages with TSS coordinates first and then peak coordinates second
write.table( tidyr::separate( tidyr::separate( with(subset(cor.mat.filt, filt.s == T), data.frame("g2" = gene.coords, "g" = gene, "p" = peak.coords, "id" = peak) ), col = "p", into = c("chr", "start", "end"), sep = "[:-]" ), col = "g2", into = c("gChr", "gStart", "gEnd"), sep = "[:-]" ), file = "atac_rna_linkages.txt", quote = F, sep = "\t", row.names = F, col.names = F )

# Load in the current genome sizes and TSS sites. 
txdb <- GenomicFeatures::makeTxDbFromGFF("mm10.refseqCurated.gtf.gz")
tss <- promoters(txdb, upstream = 0, downstream = 3)

# Get a list of all the tss regions and write them to file!
df <- data.frame(seqnames = seqnames(tss), starts = start(tss)-1, ends = end(tss)); df <- subset(df, !grepl("Un|random|alt", seqnames))
write.table(unique(df), file = "tmp", row.names = F, col.names = F, sep = "\t", quote = F)
system(paste0("sort -k1,1 -k2,2n tmp | awk '{OFS=\"\\t\"; print $0, \"tss_\" ++count}' > refSC_mm10_tss.bed"))
system("rm tmp")

# Figure out which of the TSSs overlap an ATAC-RNA linkage and only use those for the shuffle
system("intersectBed -wa -a refSC_mm10_tss.bed -b atac_rna_linkages.txt | uniq > tmp")
a <- read.table("tmp"); colnames(a) <- c("chr", "start", "end", "id")
tss <- makeGRangesFromDataFrame(a, starts.in.df.are.0based=T, keep.extra.columns=T)

# Save the "tss" object so you can use it for the shuffle 
saveRDS(tss, file = "atac_rna_cor_tss.rds")

# Now get a data.table object where each row is a tss (from above) and the entries are a list of peaks that connect with that region
system("intersectBed -wa -wb -a atac_rna_linkages.txt -b tmp | cut -f8,12 | sort | uniq > int")
lut <- as.data.table(a); setkey(lut, id)
lut$val <- list()

# The below takes a few minutes, but is still *much* faster than normal R. Hurray for data.table!
int <- read.table("int"); colnames(int) <- c("peak_id", "tss_id")
for (f in 1:nrow(int)) {
	lut[ J(int[f,"tss_id"]), val := c(lut[J(int[f,"tss_id"]),val][[1]], int[f,"peak_id"]) ] 
}

# Save the "lut" object so you don't have to regenerate it again
saveRDS(lut, file = "atac_rna_cor_lut.rds")

# Remove all the objects/files you don't need anymore
system("rm tmp int"); rm(a, int, f, df)

# Now, 
zcat /workdir/jdc397/microC/abc/immgen_cd8_rss/Predictions/EnhancerPredictionsAllPutative.tsv.gz | awk '{OFS="\t"; if ($5 != "promoter" && $28 >= 0.025) print $1, $2, $3, $4, $5, $11, $12, $16, $28}' > rss_filt_abc_links.txt

cut -f1-4 atac/summit_unif_peaks_10k.txt | intersectBed -wa -wb -a stdin -b <( sed '1d' rss_filt_abc_links.txt | cut -f1-3,6,7 ) | awk '{OFS="\t"; print $1, $9, $9, $8, $5, $6, $7, $1, $2, $3, $4}' | intersectBed -wa -wb -a stdin -b refSC_mm10_tss.bed | awk '{OFS="\t"; print $8, $9, $10, $11, int(($14 + $13 + 1) / 2), $15}' | sort -u | sort -k1,1 -k2,2n -k5,5 > rss_abc_links_int_base_shuffle.txt

module load R/4.1.3-r9
Rscript --vanilla linkageShuffle.R rss_abc_links_int_base_shuffle.txt 64

