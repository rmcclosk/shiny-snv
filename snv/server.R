library(shiny)
library(ggplot2)
library(ggvis)

shinyServer(function(input, output) {

    # selector box for patients
    make.patient.select <- function(tag, patients) {
        selectInput(tag, "Patient:", patients)
    }
    output$patient.select.vaf <- renderUI({make.patient.select("patient.vaf", unique(variants$patient))})
    output$patient.select.cna <- renderUI({make.patient.select("patient.cna", unique(segments$patient))})

    # selector boxes for samples at each time point
    # TODO: should be disabled when only one sample is available
    make.sample.select <- function (input.patient, tag) {
        m <- subset(metadata, patient==input.patient)
        do.call(tagList, as.list(by(m, m$time.point, function (s) {
            tp <- s[1, time.point]
            var <- paste0(tag, tp)
            label <- paste0("Time point ", tp, " sample:") 
            choices <- unique(as.character(s$sample))
            html <- selectInput(var, label, choices)
        })))
    }
    output$sample.select.vaf <- renderUI({make.sample.select(input$patient.vaf, "vaf.tp")})
    output$sample.select.cna <- renderUI({make.sample.select(input$patient.cna, "cna.tp")})

    # select box for chrosomes
    make.chr.select <- function (tag, selected, all) {
        html <- selectInput(tag, "Chromosome:", CHROMOSOMES, multi=T, selected=selected)
        if (all)
            gsub("(<select[^>]*)>", "\\1 disabled>", html)
        else
            html
    }
    output$chr.select.vaf <- renderUI({make.chr.select("chr.vaf", input$chr.vaf, input$all.chr.vaf)})
    output$chr.select.cna <- renderUI({make.chr.select("chr.cna", input$chr.cna, input$all.chr.cna)})

    plot.variants <- reactive({
        d <- subset(variants, patient == input$patient.vaf)

        # filter chrosomes and silent mutations
        if (!input$all.chr.vaf) d <- subset(d, chr %in% input$chr)
        if (input$hide.silent.plot) d <- subset(d, class != "Silent")
        if (nrow(d) == 0) return (add.zero.row(d))

        # remove SNVs which fall below minimum coverage depth in any sample
        min.depth <- aggregate(depth~chr+start+end+ref+alt, d, min)
        min.depth <- subset(min.depth, depth >= input$depth, select=VARIANT.COLS)
        d <- merge(d, min.depth, by=VARIANT.COLS, suffixes=c("", ".min"))

        if (input$correct.purity) d$vaf <- d$ccf

        # keep only SNVs which fall below a threshold in some sample
        if (input$maxmin < 1) {
            min.freqs <- aggregate(vaf~chr+start+end+ref+alt, d, min)
            min.freqs <- subset(min.freqs, vaf <= input$maxmin)
            d <- merge(d, min.freqs, by=VARIANT.COLS, suffixes=c("", ".min"))
        }

        samples <- sapply(unique(d$time.point), function (tp) {
            input[[paste0("vaf.tp", tp)]]
        })
        d <- d[d$sample %in% samples,]
        if (nrow(d) == 0) return (add.zero.row(d))

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
        d <- subset(segments, patient==input$patient.cna)
        samples <- sapply(unique(d$time.point), function (tp) {
            input[[paste0("cna.tp", tp)]]
        })
        d <- subset(d, sample %in% samples)
        if (nrow(d) == 0) return (add.zero.row(d))

        adj <- (1:length(samples))-floor(length(samples)/2)
        tmp <- data.frame(sample=samples, adj)

        d <- d[d$sample %in% samples,]
        d <- merge(d, tmp, by=c("sample"))
        d$copy.number <- d$copy.number + 0.15*d$adj
        d <- as.data.frame(d)
        d$abs.start <- as.numeric(d$abs.start)
        d$abs.end <- as.numeric(d$abs.end)
        d
    })

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

    # main plot
    freqPlot.vis <- reactive({
        vis <- plot.variants() %>% 
            ggvis(x=~sample, y=~vaf) %>%
            add_axis("x", title="sample") %>%
            add_axis("y", title="variant allele fraction") 

        if (input$color == "None") 
            vis <- vis %>% layer_points(key := ~key) 
        else if (input$color == "Chromosome")
            vis <- vis %>% layer_points(key := ~key, fill=~chr) %>%
                scale_nominal("fill", domain=CHROMOSOMES, range=substr(rainbow(23), 1, 7))
        else if (input$color == "Mutation type")
            vis <- vis %>% layer_points(key := ~key, fill=~class) %>%
                scale_nominal("fill", domain=unique(variants$class), range=substr(rainbow(length(unique(variants$class))), 1, 7))

        vis <- vis %>% group_by(chr, start, end, ref, alt)

        if (input$color == "None")
            vis <- vis %>% layer_paths() 
        else if (input$color == "Chromosome")
            vis <- vis %>% layer_paths(stroke=~chr) %>%
                scale_nominal("stroke", domain=CHROMOSOMES, range=substr(rainbow(23), 1, 7))
        else if (input$color == "Mutation type")
            vis <- vis %>% layer_paths(stroke=~class) %>%
                scale_nominal("stroke", domain=unique(variants$class), range=substr(rainbow(length(unique(variants$class))), 1, 7))

        vis %>% add_tooltip(tooltip, "hover")
    })

    freqPlot.vis %>% bind_shiny("freqPlot")

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

    output$table <- renderDataTable({
        overallTable()
    })

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
            geom_segment(aes(xend=abs.end, yend=copy.number, size=prevalence)) +
            theme_bw() +
            ylab("copy number") +
            xlab("chromosome") +
            theme(axis.ticks.x=element_blank()) +
            scale_x_continuous(breaks=breaks, labels=c(1:22, "X"))
    })

    # plot hierarchical clustering
    output$hclust <- renderPlot({
        n.samples <- nrow(subset(metadata, patient==input$patient & time.point > 0))
        if (n.samples > 2)
            ggdendrogram(hclusts[[input$patient]])
    })
})
