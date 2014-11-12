#!/usr/bin/env Rscript

segments <- read.table("../segments.tsv", header=T)
variants <- read.table("../snv/vaf_processed.tsv", header=T, sep="\t")

n.unique <- function (x) length(unique(x))
n.samples <- aggregate(sample~patient, segments, n.unique)
keep.patients <- subset(n.samples, sample > 1, select=c("patient"))
segments <- droplevels(merge(segments, keep.patients))
segments$chrom <- factor(segments$chrom, levels=c(1:22, "X", "Y"))

metadata <- read.table("../metadata.tsv", header=T)
segments <- merge(segments, metadata)
