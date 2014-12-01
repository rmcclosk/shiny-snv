#!/usr/bin/env Rscript

################################################################################
# Constants
################################################################################

CHROMOSOMES <- c(1:22, "X")
SNV.COLS <- c("chrom", "pos", "ref", "alt")

################################################################################
# Helper functions
################################################################################

# Find all differences in VAF between pairs of samples from different time
# points in the same patient. d is a data.frame with at least the columns
# "patient", "time.point", and a third column specified in the "col" parameter.
# Return a vector with all pairwise differences in VAF between samples of the
# same patient. For example,
#
#   patient     time.point      col
#         1              1      0.5
#         1              2      0.1
#         2              1      0.4
#         2              2      0.2
#         2              2      0.5
#
# The returned vector would contain 0.4 (for the difference between the two
# patient 1 samples), 0.2, and 0.1 (for the differences between the patient 2,
# time point 1 sample and the other two patient 2 samples). It would *not*
# include the value 0.3 (the two patient 2, time point 2 samples), because
# these samples are from the same time point.

# Note that this assumes that tumors from different time points are related, so
# that comparing them makes sense. This assumption would be violated if there
# were multiple metastases at several time points.
all.diffs <- function (d, col) {
    if (nrow(d) <= 1) return (NULL)
    first <- d[1,]
    rest <- subset(d, patient == first$patient & time.point != first$time.point)
    if (nrow(rest) > 0) {
        res <- c(rest[[col]] - first[[col]])
        names(res) <- paste0(first[["sample"]], ", ", rest[["sample"]])
        c(res, all.diffs(d[2:nrow(d),], col))
    } else {
        all.diffs(d[2:nrow(d),], col)
    }
}

# Concatenate the unique elements of x into a comma-separated string. For
# example, if x = c("a", "b", "b", "d"), then the returned string would be
# "a, b, c".
paste.unique <- function (x) { 
    paste(unique(as.character(x)), collapse=", ")
}

# Add a row of all zeros to an empty data.frame, preserving column names.
add.zero.row <- function (df) {
    names <- colnames(df)
    df <- rbind(df, rep(0, ncol(df)))
    colnames(df) <- names
    df
}

# Correct variant allele fraction based on copy number and tumor purity.
vaf.to.ccf <- function (vaf, tumor.purity, copy.number, prevalence) {
    avg.tum.copies <- copy.number*prevalence + 2*(1-prevalence)
    avg.copies <- (avg.tum.copies*tumor.purity + 2*(1-tumor.purity))
    ifelse(copy.number <= 2 | vaf <= tumor.purity / avg.copies,
           vaf * (avg.copies) / tumor.purity,
           vaf * (avg.copies) / (tumor.purity * (copy.number - 1)))
}

# find intervals which are nearest to a set of points, and return the value
# associated to each interval 
# points and intervals are assumed to be sorted
find.intervals <- function (points, intervals, value) {
    if (length(points) == 0) return (NULL)
    if (nrow(intervals) == 1) return (rep(intervals[1,value], length(points)))

    # first point is before current interval
    if (points[1] < intervals[1,1]) {
        c(intervals[1,value], find.intervals(tail(points, -1), intervals, value))

    # first point is between current and next interval
    } else if (points[1] > intervals[1,2] & points[1] < intervals[2,1]) {
        if (points[1] - intervals[1,2] <= intervals[2,1] - points[1])
            c(intervals[1,value], find.intervals(tail(points, -1), intervals, value))
        else
            c(intervals[2,value], find.intervals(tail(points, -1), tail(intervals, -1), value))

    # point is inside current interval
    } else if (points[1] >= intervals[1,1] & points[1] <= intervals[1,2]) {
        c(intervals[1,value], find.intervals(tail(points, -1), intervals, value))

    # point is at or after start of next interval
    } else {
        find.intervals(points, tail(intervals, -1), value)
    }
}

################################################################################
# Main
################################################################################

# read metadata
metadata <- read.table("../metadata.tsv", header=T)

if (file.exists("vaf_processed.tsv")) {
    d <- read.table("vaf_processed.tsv", sep="\t", header=T)
} else {
    # include only patients with at least two samples
    tumor.samples <- subset(metadata, time.point != 0)
    counts <- aggregate(sample~patient, tumor.samples, length)
    keep.patients <- unique(counts[counts$sample > 1, "patient"])
    metadata <- metadata[metadata$patient %in% keep.patients,]
    
    # read SNV information
    d <- read.table("../vaf.tsv", header=T, sep="\t", fill=T, na.strings=c("NA", ""))
    d <- d[d$chrom %in% CHROMOSOMES,]
    d$chrom <- factor(d$chrom, levels=CHROMOSOMES)
    d <- merge(d, metadata)
    d <- d[d$time.point != 0,]

    # read segments
    segments <- read.table("../segments.tsv", header=T, sep="\t")
    segments$chrom <- factor(segments$chrom, levels=CHROMOSOMES)
    
    # remove positions with depth 0 in any sample
    min.depth <- aggregate(depth~chrom+pos+ref+alt+patient, d, min)
    keep.rows <- which(min.depth$depth > 0)
    min.depth <- min.depth[keep.rows, c(SNV.COLS, "patient")]
    d <- merge(d, min.depth)
    
    # calculate variant allele frequencies
    d$vaf <- d$alt.count/d$depth
    
    # remove patients without segments
    d <- d[d$sample %in% unique(segments$sample),]

    # order by patient and sample
    d <- droplevels(d)
    d$sample <- factor(d$sample, levels=levels(segments$sample))
    d <- d[order(d$patient, d$time.point, d$sample, d$chrom, d$pos),]
    
    # give a unique key to each row
    d$key <- 1:nrow(d)

    # find copy number for each allele
    aggregate(key~sample+chrom, d, function (idx) {
        data <- d[idx,]
        by.sample <- data[1, "sample"]
        by.chrom <- data[1, "chrom"]
        intervals <- subset(segments, sample == by.sample & chrom == by.chrom,
                            select=c("start", "end", "copy.number"))
        d[idx,"copy.number"] <<- find.intervals(data$pos, intervals, "copy.number")
    })
    
    # calculate corrected VAF
    d$vaf.corrected <- vaf.to.ccf(d$vaf, d$purity/100, d$copy.number)
    conf.int <- t(mapply(function (x, n) {
        prop.test(x, n, conf.level=0.95, correct=F)$conf.int
    }, d$alt.count, d$depth))
    d$ci.lower <- vaf.to.ccf(conf.int[,1], d$purity/100, d$copy.number)
    d$ci.upper <- vaf.to.ccf(conf.int[,2], d$purity/100, d$copy.number)
    
    # find the maximum change of VAF between any pair of samples from any patient
    agg <- aggregate(key~chrom+pos+ref+alt, d, function (idx) {
        diffs <- all.diffs(d[idx,], "vaf")
        max.idx <- which.max(abs(diffs))
        c(diffs[max.idx], names(diffs)[max.idx])
    })
    agg$max.change <- agg$key[,1]
    agg$max.change.uncorrected.samples <- agg$key[,2]
    agg <- agg[,c(SNV.COLS, "max.change", "max.change.uncorrected.samples")]
    
    # same thing but with corrected VAF
    agg2 <- aggregate(key~chrom+pos+ref+alt, d, function (idx) {
        diffs <- all.diffs(d[idx,], "vaf.corrected")
        max.idx <- which.max(abs(diffs))
        c(diffs[max.idx], names(diffs)[max.idx])
    })
    agg2$max.change.corrected <- agg2$key[,1]
    agg2$max.change.corrected.samples <- agg2$key[,2]
    agg2 <- agg2[,c(SNV.COLS, "max.change.corrected", "max.change.corrected.samples")]
    
    # list all patient ID's for each SNV
    agg3 <- aggregate(patient~chrom+pos+ref+alt, d, paste.unique)
    colnames(agg3)[5] <- "patients"
    
    # merge everything
    overall <- merge(merge(agg, agg2), agg3)
    d <- merge(d, overall)
    
    write.table(d, file="vaf_processed.tsv", sep="\t", row.names=F, quote=F)
}
