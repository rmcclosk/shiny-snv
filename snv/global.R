#!/usr/bin/env Rscript

options(warn=1)
library(data.table)
library(bit64)

################################################################################
# Constants
################################################################################

CHROMOSOMES <- c(1:22, "X")
DISALLOWED.CLASSES <- c("Intron", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "IGR")
BED.COLS <- c("chr", "start", "end")
VARIANT.COLS <- c(BED.COLS, "ref", "alt")

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
    new.row <- as.list(rep(0, ncol(df)))
    names(new.row) <- colnames(df)
    rbind(df, new.row)
}

# Correct variant allele fraction based on copy number and tumor purity.
vaf.to.ccf <- function (vaf, pur, prev, minor.cn, major.cn) {
    cn <- major.cn + minor.cn
    alpha <- (cn*prev + 2*(1-prev))*pur + 2*(1-pur)

    if (minor.cn <= 1 & major.cn <= 1) {
        alpha*vaf/pur

    } else if (minor.cn == 1) {
        if (vaf >= (major.cn*prev+1)*pur/(2*alpha))
            alpha*vaf/pur - (major.cn-1)*prev
        else
            alpha*vaf/pur

    } else if (minor.cn == 0) {
        if (vaf >= (1+(major.cn-1)*prev)*pur/(2*alpha))
            alpha*vaf/pur - (major.cn-1)*prev
        else
            alpha*vaf/pur

    } else {
        if (vaf <= (1+(minor.cn-1)*prev)*pur/(2*alpha))
            alpha*vaf/pur
        else if (vaf >= (1+(major.cn+minor.cn-1)*prev)*pur/(2*alpha))
            alpha*vaf/pur - (major.cn-1)*prev
        else
            alpha*vaf/pur - (minor.cn-1)*prev
    }   
}

# rearrange a data frame to be in BED file format
# the data frame must have the columns "chr", "start", and "end"
as.bed <- function (bed.df) {
    bed.cols <- match(BED.COLS, colnames(bed.df))
    stopifnot(all(!is.na(bed.cols)))
    data.cols <- setdiff(1:ncol(bed.df), bed.cols)
    bed.df[,c(bed.cols, data.cols), with=F]
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
    bed.file.a <- write.bed(bed.a)
    cmd <- paste("bedtools", command, args)
    if (!is.null(bed.b)) {
        bed.file.b <- write.bed(bed.b)
        cmd <- paste(cmd, "-a", bed.file.a, "-b", bed.file.b)
    } else {
        cmd <- paste(cmd, "-i", bed.file.a)
    }
    cat("Running", cmd, "... ")
    output <- paste(system(cmd, intern=T), collapse="\n")
    cat("done\n")
    if (output == "") return (NULL)
    res <- fread(output)
    if (ncol(res) == 3)
        setnames(res, BED.COLS)
    else if (ncol(res) == ncol(bed.a))
        setnames(res, colnames(bed.a))
    else if (ncol(res) == ncol(bed.a) + ncol(bed.b))
        setnames(res, make.names(c(colnames(bed.a), colnames(bed.b)), unique=T))
    else if (ncol(res) == ncol(bed.a) + ncol(bed.b) + 1)
        setnames(res, make.names(c(colnames(bed.a), colnames(bed.b), "overlap"), unique=T))
    else
        res
}

# return the value in x with largest magnitude
max.abs <- function (x) {
    x[which.max(abs(x))]
}

# return the name of the value in x with largest magnitude
max.abs.name <- function (x) {
    names(x)[which.max(abs(x))]
}

################################################################################
# Main
################################################################################

# read metadata
metadata <- fread("../metadata.tsv", stringsAsFactors=T)
genome <- setNames(fread("../hg19.genome", header=F), c("chr", "chr.length"))
genome$chr.start <- cumsum(c(0, head(as.numeric(genome$chr.length), -1)))
genome <- subset(genome, chr %in% c(1:22, "X"))
metadata <- subset(metadata, time.point > 0, 
                   select=c(patient, sample, purity, time.point))

files <- c("vaf_by_patient.tsv", "vaf_overall.tsv", "cna_by_patient.tsv")
if (all(file.exists(c(files)))) {
    variants <- fread("vaf_by_patient.tsv")
    overall <- fread("vaf_overall.tsv")
    segments <- fread("cna_by_patient.tsv")
} else {
    # read SNV information
    classes <- list(character=c("chr", "rs"))
    variants <- fread("../variants.tsv", colClasses=classes)
    variants <- subset(variants, chr %in% CHROMOSOMES & ! class %in% DISALLOWED.CLASSES)

    # read segments
    segments <- fread("../segments.tsv", colClasses=list(character=c("chr")))
    segments <- subset(segments, end - start > 0)
    segments[segments$minor.cn == 1 & segments$major.cn == 1, "prevalence"] <- 1
    segments <- merge(segments, genome, by=c("chr"))
    segments$abs.start <- as.numeric(segments$chr.start + segments$start)
    segments$abs.end <- as.numeric(segments$chr.start + segments$end)

    # keep only samples with both variants and segments
    keep.samples <- intersect(unique(variants$sample), unique(segments$sample))
    metadata <- subset(metadata, sample %in% keep.samples)

    # include only patients with more than one follow-up time point
    sample.count <- aggregate(sample~patient, metadata, length)
    keep.patients <- subset(sample.count, sample > 1, select=patient)
    metadata <- merge(metadata, keep.patients, by=c("patient"))

    # merge metadata with segments and variants
    variants <- as.bed(merge(variants, metadata, by=c("sample")))
    segments <- as.bed(merge(segments, metadata, by=c("sample")))

    # find the segments which overlap each variant
    variants <- rbindlist(by(segments, segments$sample, function (by.segments) {
        # get segment overlapping each variant
        by.variants <- bedtools("sort", subset(variants, sample == by.segments[1, sample]))
        by.segments <- bedtools("sort", by.segments)
        closest <- bedtools("closest", by.variants, by.segments, args="-t first")
        closest[,c(colnames(by.variants), "copy.number", "major.cn", "minor.cn", "prevalence"), with=F]
    }, simplify=F))

    # add heterozygous segments to make prevalence up to 1
    hetero.seg <- subset(segments, major.cn != 1 | minor.cn != 1)
    hetero.seg$copy.number <- 2
    hetero.seg$major.cn <- 1
    hetero.seg$minor.cn <- 1
    hetero.seg$prevalence <- 1-hetero.seg$prevalence
    segments <- rbind(segments, hetero.seg)

    # order variants and segments
    var.order <- with(variants, order(patient, time.point, sample, chr, start))
    seg.order <- with(segments, order(patient, time.point, sample, chr, start))
    variants <- variants[var.order,]
    segments <- segments[seg.order,]

    # remove variants with depth 0 in any sample
    min.depth <- aggregate(depth~chr+start+end+ref+alt+patient, variants, min)
    min.depth <- subset(min.depth, depth > 0, select=c(VARIANT.COLS, "patient"))
    variants <- merge(variants, min.depth, by=colnames(min.depth))

    # give a unique key to each row
    variants$key <- 1:nrow(variants)
    segments$key <- 1:nrow(segments)

    # calculate variant allele frequencies
    variants$vaf <- with(variants, alt.count/depth)
    variants$ccf <- with(variants, mapply(vaf.to.ccf, vaf, purity, minor.cn, major.cn, prevalence))

    # calculate confidence intervals around variant allele frequencies
    conf.int <- t(mapply(function (x, n) {
        prop.test(x, n, conf.level=0.95, correct=F)$conf.int
    }, variants$alt.count, variants$depth))
    variants$ci.lower <- conf.int[,1]
    variants$ci.upper <- conf.int[,2]

    # find the maximum change of VAF between any pair of samples from any patient
    overall <- unique(variants[, 
        list(prot.change,
             gene,
             max.change.vaf = max.abs(all.diffs(.SD, "vaf")), 
             max.change.vaf.samples = max.abs.name(all.diffs(.SD, "vaf")), 
             max.change.ccf = max.abs(all.diffs(.SD, "ccf")),
             max.change.ccf.samples = max.abs.name(all.diffs(.SD, "ccf")),
             min.vaf = min(vaf),
             min.ccf = min(ccf),
             patients = paste.unique(patient)),
        by=list(chr,start,end,ref,alt)
    ])
    overall$start <- NULL
    overall$end <- NULL
    overall$ref <- NULL
    overall$alt <- NULL
    cat("Finished processing variants")

    write.table(overall, file="vaf_overall.tsv", sep="\t", row.names=F, quote=F)
    write.table(variants, file="vaf_by_patient.tsv", sep="\t", row.names=F, quote=F)
    write.table(segments, file="cna_by_patient.tsv", sep="\t", row.names=F, quote=F)
}

# Euclidian distance between two sets of segments
# this only takes into account segments which are overlapping between the two
# samples
sample.dist <- function (s1, s2) {
    if (s1 == s2) return (0)
    segs1 <- subset(segments, sample == s1, select=c(BED.COLS, "copy.number"))
    segs2 <- subset(segments, sample == s2, select=c(BED.COLS, "copy.number"))
    intersect <- bedtools("intersect", segs1, segs2, "-wo")
    intersect <- subset(intersect, !is.na(copy.number))
    with(intersect, sqrt(sum((abs(copy.number - copy.number.1)*overlap)^2)))
}

# hierarchical clustering for dendrograms
hclusts <- lapply(unique(segments$patient), function (by.patient) {
    samples <- subset(metadata, patient==by.patient & time.point > 0)$sample
    m <- matrix(apply(expand.grid(samples, samples), 1, function (row) {
        sample.dist(row[1], row[2])
    }), ncol=length(samples))
    rownames(m) <- samples
    colnames(m) <- samples
    hclust(as.dist(m))
})
names(hclusts) <- unique(segments$patient)
