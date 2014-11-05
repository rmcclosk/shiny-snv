#!/usr/bin/env Rscript

################################################################################
# Constants
################################################################################

CHROMOSOMES <- c(1:22, "X", "Y")
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
    if (nrow(rest) > 0)
        c(first[[col]] - rest[[col]], all.diffs(d[2:nrow(d),], col))
    else
        all.diffs(d[2:nrow(d),], col)
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
    
    # remove positions with depth 0 in any sample
    min.depth <- aggregate(depth~chrom+pos+ref+alt+patient, d, min)
    keep.rows <- which(min.depth$depth > 0)
    min.depth <- min.depth[keep.rows, c(SNV.COLS, "patient")]
    d <- merge(d, min.depth)
    
    # calculate variant allele frequencies
    d$vaf <- d$alt.count/d$depth
    
    # calculate corrected VAF by highest peak to 0.5
    peaks <- aggregate(vaf~sample, d, function (x) {
        dens <- density(x[x >= 0.1])
        dens$x[which.max(dens$y)]
    })
    d <- merge(d, peaks, by=c("sample"), suffixes=c("", ".peak"))
    d$vaf.corrected <- d$vaf/(2*d$vaf.peak)
    
    # order by patient and sample
    d <- droplevels(d)
    d <- d[order(d$patient, d$time.point, d$sample, d$chrom, d$pos),]
    
    # give a unique key to each row
    d$key <- 1:nrow(d)
    
    # find the maximum change of VAF between any pair of samples from any patient
    agg <- aggregate(key~chrom+pos+ref+alt, d, function (idx) {
        diffs <- all.diffs(d[idx,], "vaf")
        diffs[which.max(abs(diffs))]
    })
    colnames(agg)[5] <- "max.change"
    
    # same thing but with corrected VAF
    agg2 <- aggregate(key~chrom+pos+ref+alt, d, function (idx) {
        diffs <- all.diffs(d[idx,], "vaf.corrected")
        diffs[which.max(abs(diffs))]
    })
    colnames(agg2)[5] <- "max.change.corrected"
    
    # list all patient ID's for each SNV
    agg3 <- aggregate(patient~chrom+pos+ref+alt, d, paste.unique)
    colnames(agg3)[5] <- "patients"
    
    # merge everything
    overall <- merge(merge(agg, agg2), agg3)
    d <- merge(d, overall)
    
    write.table(d, file="vaf_processed.tsv", sep="\t", row.names=F, quote=F)
}
