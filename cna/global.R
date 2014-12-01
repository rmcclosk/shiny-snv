#!/usr/bin/env Rscript

segments <- read.table("../segments.tsv", header=T)
variants <- read.table("../snv/vaf_processed.tsv", header=T, sep="\t")
metadata <- read.table("../metadata.tsv", header=T)
segments <- merge(segments, metadata)

n.unique <- function (x) length(unique(x))
n.samples <- aggregate(sample~patient, segments, n.unique)
keep.patients <- subset(n.samples, sample > 1, select=c("patient"))
segments <- droplevels(merge(segments, keep.patients))
segments$chrom <- factor(segments$chrom, levels=c(1:22, "X", "Y"))

chr.ends <- aggregate(end~chrom, segments, max)
chr.ends$chr.start <- cumsum(c(0, tail(as.numeric(chr.ends$end), -1)))
chr.ends <- chr.ends[,c("chrom", "chr.start")]

segments <- merge(segments, chr.ends)
segments$adj.end <- as.numeric(segments$chr.start + segments$end)
segments$adj.start <- as.numeric(segments$chr.start + segments$start)

variants <- merge(variants, chr.ends)
variants$adj.pos <- as.numeric(variants$chr.start + variants$pos)

sample.dist <- function (s1, s2) {
    segs1 <- subset(segments, sample == s1, select=c(adj.start, adj.end, copy.number))
    colnames(segs1)[colnames(segs1) == "copy.number"] <- "copy.number.1"
    segs2 <- subset(segments, sample == s2, select=c(adj.start, adj.end, copy.number))
    colnames(segs2)[colnames(segs2) == "copy.number"] <- "copy.number.2"
    merged.segs <- merge(segs1, segs2)
    len <- merged.segs$adj.end-merged.segs$adj.start
    diff <- abs(merged.segs$copy.number.1 - merged.segs$copy.number.2)
    sqrt(sum((diff*len)^2))
}

# hierarchical clustering for dendrograms
hclusts <- lapply(levels(segments$patient), function (by.patient) {
    samples <- subset(metadata, patient==by.patient & time.point > 0)$sample
    m <- matrix(apply(expand.grid(samples, samples), 1, function (row) {
        sample.dist(row[1], row[2])
    }), ncol=length(samples))
    rownames(m) <- samples
    colnames(m) <- samples
    hclust(as.dist(m))
})
names(hclusts) <- levels(segments$patient)
