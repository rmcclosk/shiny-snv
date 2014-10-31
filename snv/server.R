library(shiny)
library(ggvis)

################################################################################
# Constants
################################################################################

CHROMOSOMES <- c(1:22, "X", "Y")

# columns which uniquely determine SNVs
SNV.COLS <- c("chrom", "pos", "ref", "alt")


################################################################################
# Helper functions
################################################################################

# Find all differences in VAF between pairs of samples from different time
# points in the same patient. d is a data.frame with at least the columns
# "patient", "time.point", and "vaf". Return a vector with all pairwise
# differences in VAF between samples of the same patient. For example,
#
#   patient     time.point      vaf
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
all.diffs <- function (d) {
    if (nrow(d) <= 1) return (NULL)
    first <- d[1,]
    rest <- subset(d, patient == first$patient & time.point != first$time.point)
    if (nrow(rest) > 0)
        c(first$vaf - rest$vaf, all.diffs(d[2:nrow(d),]))
    else
        all.diffs(d[2:nrow(d),])
}

# Concatenate the unique elements of x into a comma-separated string. For
# example, if x = c("a", "b", "b", "d"), then the returned string would be
# "a, b, c".
unique.list <- function (x) paste(unique(as.character(x)), collapse=", ")

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

# include only patients with at least two tumor samples
tumor.samples <- subset(metadata, time.point != 0)
counts <- aggregate(sample~patient, tumor.samples, length)
keep.patients <- unique(counts[counts$sample > 1, "patient"])
metadata <- metadata[metadata$patient %in% keep.patients,]

# read SNV information
d <- read.table("../snv.tsv", header=T, sep="\t", fill=T, na.strings=c("NA", ""))
d <- d[d$chrom %in% CHROMOSOMES,]
d$chrom <- factor(d$chrom, levels=CHROMOSOMES)
d <- merge(d, metadata)

# remove positions with depth 0 in any sample
min.depth <- aggregate(depth~chrom+pos+ref+alt+patient, d, min)
keep.rows <- which(min.depth$depth > 0)
min.depth <- min.depth[keep.rows, c("chrom", "pos", "ref", "alt", "patient")]
d <- merge(d, min.depth)

# calculate variant allele frequencies
d$vaf <- d$alt.count/d$depth

# order by patient and sample
d <- droplevels(d)
d <- d[order(d$patient, d$time.point, d$sample, d$chrom, d$pos),]

# give a unique key to each row
d$key <- 1:nrow(d)

# find the maximum change of VAF between any pair of samples from any patient
agg <- aggregate(key~chrom+pos+ref+alt, d, function (idx) {
    diffs <- all.diffs(d[idx,])
    diffs[which.max(abs(diffs))]
})
colnames(agg)[5] <- "max.change"

# list all patient ID's for each SNV
agg2 <- aggregate(patient~chrom+pos+ref+alt, d, unique.list)
colnames(agg2)[5] <- "patients"
overall <- merge(agg, agg2)

# combine everything into a single table
keep.cols <- c("chrom", "pos", "ref", "alt", "prot.change", "gene", "class", 
               "max.change", "patients")
overall <- merge(d, overall)[,keep.cols]
overall <- overall[!duplicated(overall),]

shinyServer(function(input, output) {

    # selector box for patients
    output$patient.select <- renderUI({
        selectInput("patient", "Patient:", levels(d$patient))
    })

    # selector boxes for samples at each time point
    # TODO: should be disabled when only one sample is available
    output$sample.select <- renderUI({
        m <- subset(metadata, patient==input$patient & time.point != 0)
        do.call(tagList, as.list(by(m, m$time.point, function (s) {
            tp <- s[1, "time.point"]
            var <- paste0("tp", tp)
            label <- paste0("Time point ", tp, " sample:") 
            choices <- unique(as.character(s$sample))
            html <- selectInput(var, label, choices)
        })))
    })

    # select box for chromosomes
    output$chrom.select <- renderUI({
        html <- selectInput("chrom", "Chromosome:", levels(d$chrom), multi=T, selected=input$chrom)
        if (input$all.chrom) 
            gsub("(<select[^>]*)>", "\\1 disabled>", html)
        else
            html
    })

    plot.data <- reactive({
        pd <- subset(d, patient == input$patient & time.point != 0)

        # filter chromosomes
        if (!input$all.chrom) pd <- subset(pd, chrom %in% input$chrom)

        # filter silent mutations
        if (input$hide.silent) pd <- subset(pd, class != "Silent")

        if (nrow(pd) == 0) return (add.zero.row(pd))

        # remove SNVs which fall below minimum coverage depth in any sample
        min.depth <- aggregate(depth~chrom+pos+ref+alt, pd, min)
        min.depth <- subset(min.depth, depth >= input$depth, select=SNV.COLS)
        pd <- merge(pd, min.depth, by=SNV.COLS, suffixes=c("", ".min"))

        # keep only SNVs which fall below a threshold in some sample
        min.freqs <- aggregate(vaf~chrom+pos+ref+alt, pd, min)
        min.freqs <- subset(min.freqs, vaf <= input$maxmin)
        pd <- merge(pd, min.freqs, by=SNV.COLS, suffixes=c("", ".min"))

        samples <- sapply(unique(pd$time.point), function (tp) {
            input[[paste0("tp", tp)]]
        })
        pd <- pd[pd$sample %in% samples,]

        if (nrow(pd) == 0) return (add.zero.row(pd))

	    if (input$order == "Highest fraction") {
	        aggfun <- sum
	    } else if (input$order == "Most change") {
	        aggfun <- function (x) max(abs(diff(x)))
	    }
	    agg <- aggregate(vaf~chrom+pos+ref+alt, pd, aggfun)
	    agg <- agg[order(-agg$vaf),]
	    agg <- head(agg, input$n)
        pd <- merge(pd, agg, by=SNV.COLS, suffixes=c("", ".agg"))

        pd <- droplevels(pd[order(pd$time.point),])
        if (input$correct.purity) pd$vaf <- pd$vaf/(pd$purity/100)
        pd
    })

    tooltip <- function (x) {
        pd <- isolate(plot.data())
        point <- pd[pd$key == x$key,]
        if (nrow(point) == 0) return ("")
        html <- paste0("<b>Chromosome: </b>", point$chrom, "</b><br />")
        html <- paste0(html, "<b>Position: </b>", point$pos, "</b><br />")
        html <- paste0(html, "<b>Reference: </b>", point$ref, "</b><br />")
        html <- paste0(html, "<b>Variant: </b>", point$alt, "</b><br />")
        html <- paste0(html, "<b>Gene: </b>", point$gene, "</b><br />")
        if (!is.na(point$prot.change))
            html <- paste0(html, "<b>Protein change: </b>", point$prot.change, "</b><br />")
        if (!is.na(point$rs))
            html <- paste0(html, "<b>dbSNP ID: </b>", point$rs, "</b><br />")
        if (!is.na(point$cosmic))
            html <- paste0(html, "<b>COSMIC ID: </b>", point$cosmic, "</b><br />")
        if (!is.na(point$esp))
            html <- paste0(html, "<b>ESP ID: </b>", point$esp, "</b><br />")
        html
    }

    # main plot
    freqPlot.vis <- reactive({
        vis <- plot.data %>% 
            ggvis(x=~factor(sample), y=~vaf) %>%
            add_axis("x", title="sample") %>%
            add_axis("y", title="variant allele fraction") 

        if (input$color == "None") 
            vis <- vis %>% layer_points(key := ~key) 
        else if (input$color == "Chromosome")
            vis <- vis %>% layer_points(key := ~key, fill=~chrom)
        else if (input$color == "Mutation type")
            vis <- vis %>% layer_points(key := ~key, fill=~class)

        vis <- vis %>% group_by(chrom, pos, ref, alt)

        if (input$color == "None")
            vis <- vis %>% layer_paths() 
        else if (input$color == "Chromosome")
            vis <- vis %>% layer_paths(stroke=~chrom)
        else if (input$color == "Mutation type")
            vis <- vis %>% layer_paths(stroke=~class)

        vis %>% add_tooltip(tooltip, "hover")
    })

    freqPlot.vis %>% bind_shiny("freqPlot")

    output$freqTable <- renderDataTable({
        overall
    })
})
