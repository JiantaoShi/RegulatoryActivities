library("rGREAT")
library("dplyr")

df = read.table(file = "GeneTSS_Hg19.tsv", row.names = NULL, header = TRUE, sep = "\t")
tss = GRanges(seqnames = df[, 2], ranges = IRanges(df[, 3], df[, 3]), strand = df[, 4], gene_id = df[,1], gene_name = df[,5], interval = df[,6])

# remove duplicate TSS
smx = summary(factor(tss$gene_id), maxsum = length(tss))
rmid = names(smx)[smx > 1] # 1 duplicates
t0 = tss[!(tss$gene_id %in% rmid)]
t1 = tss[tss$gene_id %in% rmid]
tss = sort(c(t0, t1[1]))
rLocus = extendTSS(tss, genome = "hg19", gene_id_type = "ENSEMBL")

# promoter
newTSS = GRanges(seqnames = seqnames(rLocus), ranges = IRanges(rLocus$tss_position, rLocus$tss_position), 
	strand = rLocus$tss_strand, gene_id = rLocus$gene_id)
pLocus = promoters(newTSS, upstream=5000, downstream=1000)

# out
names(tss) = tss$gene_id
out = data.frame(gene_id = pLocus$gene_id,
				CHROM = seqnames(pLocus),
				POS = rLocus$tss_position,
				STRAND = strand(pLocus),
				gene_name = tss[pLocus$gene_id]$gene_name,
				gene_locus = tss[pLocus$gene_id]$interval,
				reg_locus = paste0(seqnames(rLocus), ":", start(rLocus), "-", end(rLocus)),
				promoter = paste0(seqnames(pLocus), ":", start(pLocus), "-", end(pLocus))
				)

write.table(out, file = "GeneAnnotation_Hg19.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE)
