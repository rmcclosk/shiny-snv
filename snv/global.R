#!/usr/bin/env Rscript

################################################################################
# Constants
################################################################################

CHROMOSOMES <- c(1:22, "X")
DISALLOWED.CLASSES <- c("Intron", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "IGR")
BED.COLS <- c("chr", "start", "end")
VARIANT.COLS <- c(BED.COLS, "start", "end")

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
    if (copy.number <= 2)
        return(vaf*avg.copies/tumor.purity)

    ifelse(copy.number <= 2 | vaf <= tumor.purity / avg.copies,
           vaf * (avg.copies) / tumor.purity,
           vaf * (avg.copies) / (tumor.purity * (copy.number - 1)))
}

# rearrange a data frame to be in BED file format
# the data frame must have the columns "chr", "start", and "end"
as.bed <- function (bed.df) {
    bed.cols <- match(BED.COLS, colnames(bed.df))
    stopifnot(all(!is.na(bed.cols)))
    data.cols <- setdiff(1:ncol(bed.df), bed.cols)
    bed.df[,c(bed.cols, data.cols)]
}

# write a data frame to a BED file
# data frame must have the columns "chr", "start", and "end"
# return the file name written to
write.bed <- function (bed.df, bed.file=tempfile()) {
    stopifnot(colnames(bed.df)[1:length(BED.COLS)] == BED.COLS)
    write.table(bed.df, bed.file, col.names=F, row.name=F, sep="\t", quote=F)
    bed.file
}

# run a bedtools command 
# bed files should both have chr, start, end as their first three columns
bedtools <- function (command, bed.a, bed.b=NULL, args="") {
    bed.file.a <- write.bed(bed.a, "a.bed")
    cmd <- paste("bedtools", command, args)
    if (!is.null(bed.b)) {
        bed.file.b <- write.bed(bed.b, "b.bed")
        cmd <- paste(cmd, "-a", bed.file.a, "-b", bed.file.b)
    } else {
        cmd <- paste(cmd, "-i", bed.file.a)
    }
    cat("Running", cmd, "... ")
    output <- system(cmd, intern=T)
    cat("done\n")
    if (length(output) == 0) return (NULL)
    read.table(textConnection(output), sep="\t")
}

################################################################################
# Main
################################################################################

# read metadata
metadata <- read.table("../metadata.tsv", header=T)

if (file.exists("vaf_processed.tsv")) {
    d <- read.table("vaf_processed.tsv", sep="\t", header=T)
} else {
    # read SNV information
    variants <- read.table("../variants.tsv", header=T, sep="\t")
    variants <- subset(variants, chr %in% CHROMOSOMES & ! class %in% DISALLOWED.CLASSES)
    variants$chr <- factor(variants$chr, levels=CHROMOSOMES)

    # read segments
    segments <- read.table("../segments.tsv", header=T, sep="\t")
    segments$chr <- factor(segments$chr, levels=CHROMOSOMES)
    segments <- subset(segments, end - start > 0)
    segments[is.na(segments$prevalence), "prevalence"] <- 1

    # include only samples with both segments and variants
    keep.samples <- intersect(levels(segments$sample), levels(variants$sample))
    metadata <- subset(metadata, sample %in% keep.samples)

    # include only patients with more than one follow-up time point
    counts <- aggregate(sample~patient, metadata, length)
    keep.patients <- counts[counts$sample > 1, "patient"]
    metadata <- droplevels(metadata[metadata$patient %in% keep.patients,])

    # drop samples failing any previous criteria
    keep.samples <- metadata$sample
    print(keep.samples)
    variants <- droplevels(as.bed(subset(variants, sample %in% keep.samples)))
    levels(variants$sample) <- keep.samples
    segments <- droplevels(as.bed(subset(segments, sample %in% keep.samples)))
    levels(segments$sample) <- keep.samples

    # find the segments which overlap each variant
    variants <- do.call(rbind, by(segments, segments$sample, function (by.segments) {
        # fill in gaps in the segments
        compl <- bedtools("complement", by.segments, args="-g ../hg19.genome")
        colnames(compl) <- BED.COLS
        compl$sample <- by.segments[1, "sample"]
        compl$copy.number <- 2
        compl$prevalence <- 1
        by.segments <- rbind(by.segments, compl)

        by.variants <- subset(variants, sample == by.segments[1, "sample"])
        intersect <- bedtools("intersect", by.variants, by.segments, "-loj")
        old.names <- c(colnames(by.variants), colnames(by.segments))
        new.names <- make.names(old.names, unique=T)
        colnames(intersect) <- new.names
        intersect[,c(colnames(by.variants), "copy.number", "prevalence")]
    }, simplify=F))

    # remove variants with depth 0 in any sample
    min.depth <- aggregate(depth~chr+start+end+ref+alt+patient, variants, min)
    min.depth <- subset(min.depth, depth > 0, select=c(VARIANT.COLS, "patient"))
    variants <- merge(variants, min.depth)

    # order variants and segments
    var.order <- with(variants, order(patient, time.point, sample, chr, start))
    seg.order <- with(segments, order(patient, time.point, sample, chr, start))
    variants <- variants[var.order,]
    segments <- segments[seg.order,]

    # give a unique key to each row
    variants$key <- 1:nrow(variants)
    segments$key <- 1:nrow(segments)

    # calculate variant allele frequencies
    variants$vaf <- variants$alt.count/variants$depth
    
    # calculate corrected VAF
    variants$vaf.corrected <- vaf.to.ccf(variants$vaf, variants$purity, variants$copy.number, variants$prevalence)
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
    agg <- agg[,c(VARIANT.COLS, "max.change", "max.change.uncorrected.samples")]
    
    # same thing but with corrected VAF
    agg2 <- aggregate(key~chrom+pos+ref+alt, d, function (idx) {
        diffs <- all.diffs(d[idx,], "vaf.corrected")
        max.idx <- which.max(abs(diffs))
        c(diffs[max.idx], names(diffs)[max.idx])
    })
    agg2$max.change.corrected <- agg2$key[,1]
    agg2$max.change.corrected.samples <- agg2$key[,2]
    agg2 <- agg2[,c(VARIANT.COLS, "max.change.corrected", "max.change.corrected.samples")]
    
    # list all patient ID's for each SNV
    agg3 <- aggregate(patient~chrom+pos+ref+alt, d, paste.unique)
    colnames(agg3)[5] <- "patients"
    
    # merge everything
    overall <- merge(merge(agg, agg2), agg3)
    d <- merge(d, overall)
    
    write.table(d, file="vaf_processed.tsv", sep="\t", row.names=F, quote=F)
}
