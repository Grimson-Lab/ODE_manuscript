

suppressPackageStartupMessages({ 
	library(regioneR)
})


txdb <- GenomicFeatures::makeTxDbFromGFF("../mm10.refseqCurated.gtf")
tss <- promoters(txdb, upstream = 3, downstream = 0)
write.table( unique(as.data.frame(tss)[,c(1:3)]), file = "refSeqCurated_tss.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Now have a version that is promoter-proximal, which is close enough to be considered part of the likely regulatory module
tss <- promoters(txdb, upstream = 2000, downstream = 2000)
write.table( unique(as.data.frame(tss)[,c(1:3)]), file = "refSeqCurated_tssProximal.bed", row.names = F, col.names = F, quote = F, sep = "\t")


# Figure out how many peaks in total overlap TSS sites from the above file (out of 332869 total peaks)
system(" sort -k1,1 -k2,2n refSeqCurated_tss.bed | intersectBed -sorted -u -wa -a summit_unif_peaks_10k.bed -b stdin | wc -l ")
# 16497

# How many of the peaks are promoter promixal?
system(" sort -k1,1 -k2,2n refSeqCurated_tssProximal.bed | intersectBed -sorted -u -wa -a summit_unif_peaks_10k.bed -b stdin | wc -l ")
# 33894

# How many of the accessible TEs directly overlap a TSS site (out of 123144 total accessible TEs)?
system(" sort -k1,1 -k2,2n refSeqCurated_tss.bed | intersectBed -sorted -u -wa -a accessible.tes.atac.txt -b stdin | wc -l ")
# 302

# How many of the accessible TEs are promoter promixal?
system(" sort -k1,1 -k2,2n refSeqCurated_tssProximal.bed | intersectBed -sorted -u -wa -a accessible.tes.atac.txt -b stdin | wc -l ")
# 5969


# How many of the peaks associated with accessible TEs intersect with TSS sites (out of 113340 total peaks)?
system(" sort -k1,1 -k2,2n refSeqCurated_tss.bed | intersectBed -sorted -u -wa -a peak.te.int -b stdin | cut -f1-4 | uniq | wc -l ")
# 1705





 