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

# find all differences in VAF between pairs of samples from different time
# points in the same patient
all.diffs <- function (d) {
    if (nrow(d) == 1) return (NULL)
    first <- d[1,]
    rest <- subset(d, patient == first$patient & time.point != first$time.point)
    if (nrow(rest) > 0)
        c(first$vaf - rest$vaf, all.diffs(d[2:nrow(d),]))
    else
        all.diffs(d[2:nrow(d),])
}

idx <- 1:nrow(d)
agg <- aggregate(idx~chrom+pos+ref+alt, d, function (idx) {
    if(length(idx) < 2) return(NA)
    diffs <- all.diffs(d[idx,])
    diffs[which.max(abs(diffs))]
})
colnames(agg)[5] <- "max.change"

agg2 <- aggregate(patient~chrom+pos+ref+alt, d, length)
colnames(agg2)[5] <- "n.patients"
agg <- merge(agg, agg2)
head(agg, 10)

# todo: investigate NA's
agg <- agg[!is.na(agg$max.change),]
keep.cols <- c("chrom", "pos", "ref", "alt", "prot.change", "gene", "class", 
               "max.change", "n.patients")
agg <- merge(d, agg)[,keep.cols]
agg <- agg[!duplicated(agg),]
