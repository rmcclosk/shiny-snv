library(shiny)
library(ggvis)

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

        if (input$correct.purity) pd$vaf <- pd$vaf.corrected

        # keep only SNVs which fall below a threshold in some sample
        if (input$maxmin < 1) {
            min.freqs <- aggregate(vaf~chrom+pos+ref+alt, pd, min)
            min.freqs <- subset(min.freqs, vaf <= input$maxmin)
            pd <- merge(pd, min.freqs, by=SNV.COLS, suffixes=c("", ".min"))
        }

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
        pd
    })

    tooltip <- function (x) {
        pd <- isolate(plot.data())
        if (is.null(x$key)) {
            point <- subset(pd, chrom==x$chrom & pos==x$pos & ref==x$ref & alt==x$alt)[1,]
        } else {
            point <- pd[pd$key == x$key,]
        }
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

    overallTable <- reactive({
        keep.cols <- c("chrom", "pos", "ref", "alt", "prot.change", "gene",
                       "class", "max.change", "max.change.corrected",
                       "patients")
        if (input$hide.silent.table)
            ft <- subset(d, class != "Silent", select=keep.cols)
        else
            ft <- subset(d, select=keep.cols)
        ft[!duplicated(ft),]
    })

    output$freqTable <- renderDataTable({
        overallTable()
    })

    output$downloadData <- downloadHandler(
        filename = function() { "data.csv" },
        content = function(file) { write.csv(overallTable(), file) }
    )
})
