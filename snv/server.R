library(shiny)
library(ggplot2)
library(ggvis)

# make a select input, possibly disabled
make.select <- function(var, lbl, choices, multi=F, sel=NULL, disabled=F) {
    html <- selectInput(var, lbl, choices, selected=sel, multiple=multi)
    if (!is.null(disabled) && disabled)
        gsub("(<select[^>]*)>", "\\1 disabled>", html)
    else
        html
}

# make selector boxes for samples at each time point for a patient
# input.patient: patient to make selectons for
# tag: prefix for input ID's for selectors (eg. selectors will assign to
#      variables tag1, tag2 for time points 1 and 2)
make.sample.select <- function (input.patient, tag) {
    m <- subset(metadata, patient == input.patient)
    do.call(tagList, by(m, m$time.point, function (s) {
        tp <- s[1, time.point]
        var <- paste0(tag, tp)
        lbl <- paste0("Time point ", tp, " sample:") 
        choices <- unique(s$sample)
        html <- make.select(var, lbl, choices, disabled=length(choices)==1)
    }, simplify=FALSE))
}

# constants
CLASSES <- unique(variants$class)

# colors for plots
chr.cols <- substr(rainbow(23), 1, 7)
class.cols <- substr(rainbow(length(unique(CLASSES))), 1, 7)

shinyServer(function(input, output) {

    # sample selectors
    output$sample.select.vaf <- renderUI({
        make.sample.select(input$patient.vaf, "vaf.tp")
    })
    output$sample.select.cna <- renderUI({
        make.sample.select(input$patient.cna, "cna.tp")
    })

    # chromosome and gene selectors
    output$chr.select <- renderUI({
        make.select("chr", "Chromosome:", CHROMOSOMES, multi=TRUE,
                    sel=input$chr, disabled=input$all.chr)
    })
    output$gene.select <- renderUI({
        genes <- unique(subset(variants, patient == input$patient.vaf)$gene)
        make.select("gene", "Gene:", genes, multi=TRUE,
                    sel=input$gene, disabled=input$all.gene)
    })

    plot.variants <- reactive({
        d <- subset(variants, patient == input$patient.vaf)

        # only variants from selected samples
        samples <- sapply(unique(d$time.point), function (tp) {
            input[[paste0("vaf.tp", tp)]]
        })
        d <- d[d$sample %in% samples,]

        # filter chrosomes and silent mutations
        if (!input$all.chr) d <- subset(d, chr %in% input$chr)
        if (input$hide.silent.plot) d <- subset(d, class != "Silent")

        # filter genes
        if (!input$all.gene && input$color != "Selected genes")
            d <- subset(d, gene %in% input$gene)
        if (nrow(d) == 0) return (add.zero.row(d))

        # remove SNVs which fall below minimum coverage depth in any sample
        min.depth <- aggregate(depth~chr+start+end+ref+alt, d, min)
        min.depth <- subset(min.depth, depth >= input$depth, select=VARIANT.COLS)
        d <- merge(d, min.depth, by=VARIANT.COLS, suffixes=c("", ".min"))

        # correct VAF to CCF
        if (input$correct.purity) d$vaf <- d$ccf

        # keep only SNVs which fall below a threshold in some sample
        if (input$maxmin < 1) {
            min.freqs <- aggregate(vaf~chr+start+end+ref+alt, d, min)
            min.freqs <- subset(min.freqs, vaf <= input$maxmin)
            d <- merge(d, min.freqs, by=VARIANT.COLS, suffixes=c("", ".min"))
        }
        if (nrow(d) == 0) return (add.zero.row(d))

        # order variants and take the top n
	    if (input$order == "Highest fraction") {
	        aggfun <- sum
	    } else if (input$order == "Most change") {
	        aggfun <- function (x) max(abs(diff(x)))
	    }
	    agg <- aggregate(vaf~chr+start+end+ref+alt, d, aggfun)
	    agg <- agg[order(-agg$vaf),]
	    agg <- head(agg, input$n)
        d <- merge(d, agg, by=VARIANT.COLS, suffixes=c("", ".agg"))

        as.data.frame(d[order(d$time.point, d$chr, d$start),])
    })

    # get segments to plot
    plot.segments <- reactive({

        # get segments from selected samples only
        d <- subset(segments, patient==input$patient.cna)
        samples <- sapply(unique(d$time.point), function (tp) {
            input[[paste0("cna.tp", tp)]]
        })
        d <- subset(d, sample %in% samples)
        if (nrow(d) == 0) return (add.zero.row(d))

        # y-axis adjustment (so that the segments don't overlap on the plot)
        adj <- (1:length(samples))-floor(length(samples)/2)
        tmp <- data.frame(sample=samples, adj)
        d <- merge(d, tmp, by=c("sample"))
        d$copy.number <- d$copy.number + 0.15*d$adj

        # annoying adjustments, because ggplot doesn't like data.tables
        d <- as.data.frame(d)
        d$abs.start <- as.numeric(d$abs.start)
        d$abs.end <- as.numeric(d$abs.end)
        d
    })

    # mouseover for VAF plot
    tooltip <- function (x) {
        d <- isolate(plot.variants())
        if (is.null(x$key)) {
            point <- subset(d, chr==x$chr & start==x$start & end==x$end & ref==x$ref & alt==x$alt)[1,]
        } else {
            point <- d[d$key == x$key,]
        }
        if (nrow(point) == 0) return ("")
        html <- paste0("<b>Chromosome: </b>", point$chr, "</b><br />")
        html <- paste0(html, "<b>Gene: </b>", point$gene, "</b><br />")
        html <- paste0(html, "<b>Base change: </b>", point$ref, ">", point$alt, "</b><br />")
        if (!is.na(point$prot.change) & point$prot.change != "")
            html <- paste0(html, "<b>Protein change: </b>", point$prot.change, "</b><br />")
        if (!is.na(point$rs) & point$rs != "")
            html <- paste0(html, "<b>dbSNP ID: </b>", point$rs, "</b><br />")
        if (!is.na(point$cosmic) & point$cosmic != "")
            html <- paste0(html, "<b>COSMIC ID: </b>", point$cosmic, "</b><br />")
        if (!is.na(point$esp) & point$esp != "")
            html <- paste0(html, "<b>ESP ID: </b>", point$esp, "</b><br />")
        html
    }

    # VAF plot
    freqPlot.vis <- reactive({
        vis <- plot.variants() %>% 
            ggvis(x=~sample, y=~vaf) %>%
            add_axis("x", title="sample") %>%
            add_axis("y", title="variant allele fraction") 

        if (input$color == "None") 
            vis <- vis %>% layer_points(key := ~key) 
        else if (input$color == "Chromosome")
            vis <- vis %>% layer_points(key := ~key, fill=~chr) %>%
                scale_nominal("fill", domain=CHROMOSOMES, range=chr.cols)
        else if (input$color == "Mutation type")
            vis <- vis %>% layer_points(key := ~key, fill=~class) %>%
                scale_nominal("fill", domain=CLASSES, range=class.cols)

        vis <- vis %>% group_by(chr, start, end, ref, alt)

        if (input$color == "None")
            vis <- vis %>% layer_paths() 
        else if (input$color == "Chromosome")
            vis <- vis %>% layer_paths(stroke=~chr) %>%
                scale_nominal("stroke", domain=CHROMOSOMES, range=chr.cols)
        else if (input$color == "Mutation type")
            vis <- vis %>% layer_paths(stroke=~class) %>%
                scale_nominal("stroke", domain=CLASSES, range=class.cols)

        vis %>% add_tooltip(tooltip, "hover")
    })
    freqPlot.vis %>% bind_shiny("freqPlot")

    # table for VAF overall
    overallTable <- reactive({
        tbl <- overall
        if (input$hide.silent.table)
            tbl <- subset(tbl, prot.change != "")

        tbl$max.change.vaf <- round(tbl$max.change.vaf, 3)
        tbl$max.change.ccf <- round(tbl$max.change.ccf, 3)
        tbl$min.ccf <- round(tbl$min.ccf, 3)
        tbl$min.vaf <- round(tbl$min.vaf, 3)
        tbl
    })
    output$table <- renderDataTable({ overallTable() })

    # download VAF overall table
    output$downloadData <- downloadHandler(
        filename = function() { "data.csv" },
        content = function(file) { write.csv(overallTable(), file) }
    )

    # plot segments
    output$segsPlot <- renderPlot({
        d <- plot.segments()
        starts <- genome$chr.start
        breaks <- starts + genome$chr.length/2
        ggplot(d, aes(x=abs.start, y=copy.number, color=sample)) +
            scale_x_continuous(limits=c(0, sum(as.numeric(genome$chr.length))), expand=c(0, 0), breaks=breaks, labels=c(1:22, "X")) +
            geom_segment(aes(xend=abs.end, yend=copy.number, size=prevalence)) +
            theme_bw() +
            ylab("copy number") +
            xlab("chromosome") +
            theme(axis.ticks.x=element_blank()) +
            geom_vline(data=genome, aes(xintercept=chr.start), linetype="dashed", color="grey")
    })

    # plot hierarchical clustering
    output$hclust <- renderPlot({
        n.samples <- nrow(subset(metadata, patient==input$patient & time.point > 0))
        if (n.samples > 2)
            ggdendrogram(hclusts[[input$patient]])
    })
})
