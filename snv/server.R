library(shiny)
library(ggvis)

shinyServer(function(input, output) {

    # selector box for patients
    output$patient.select <- renderUI({
        selectInput("patient", "Patient:", unique(variants$patient))
    })

    # selector boxes for samples at each time point
    # TODO: should be disabled when only one sample is available
    output$sample.select <- renderUI({
        m <- subset(metadata, patient==input$patient)
        do.call(tagList, as.list(by(m, m$time.point, function (s) {
            tp <- s[1, time.point]
            var <- paste0("tp", tp)
            label <- paste0("Time point ", tp, " sample:") 
            choices <- unique(as.character(s$sample))
            html <- selectInput(var, label, choices)
        })))
    })

    # select box for chrosomes
    output$chr.select <- renderUI({
        html <- selectInput("chr", "Chromosome:", CHROMOSOMES, multi=T, selected=input$chr)
        if (input$all.chr) 
            gsub("(<select[^>]*)>", "\\1 disabled>", html)
        else
            html
    })

    plot.variants <- reactive({
        d <- subset(variants, patient == input$patient)

        # filter chrosomes and silent mutations
        if (!input$all.chr) d <- subset(d, chr %in% input$chr)
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
            input[[paste0("tp", tp)]]
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
})
