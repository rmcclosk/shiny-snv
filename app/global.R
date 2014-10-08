#!/usr/bin/env Rscript

sample.info <- read.table("metadata.csv", header=T)
sample.info$date <- as.Date(sample.info$date, "%m/%d/%Y")

counts <- aggregate(tumor.sample~patient, sample.info, length)
sample.info <- merge(sample.info, counts, by=c("patient"), 
                     suffixes=c("", ".count"))
sample.info <- sample.info[sample.info$tumor.sample.count > 1,]

d <- read.table("freqs.dat", header=T, sep="\t", fill=T, na.strings=c("NA", ""))
d <- d[d$chrom %in% c(1:22, "X", "Y"),]
d$chrom <- factor(d$chrom, levels=c(1:22, "X", "Y"))
d <- d[d$depth > 0,]
d$vaf <- d$alt.count/d$depth

d <- merge(d, sample.info, by.x=c("patient", "sample"),
                           by.y=c("patient", "tumor.sample"))

# keep patients and samples in the same order
d <- droplevels(d)
d <- d[order(d$patient, d$date, d$sample, d$chrom, d$pos),]
d$sample <- factor(d$sample, levels=unique(d$sample))

# give a unique key so we can make a tooltip
d$key <- 1:nrow(d)
