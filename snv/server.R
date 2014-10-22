library(shiny)
library(ggvis)

add.zero.row <- function (df) {
    names <- colnames(df)
    df <- rbind(df, rep(0, ncol(df)))
    colnames(df) <- names
    df
}

shinyServer(function(input, output) {

    plot.data <- reactive({
        pd <- subset(d, patient == input$patient)
        min.depth <- aggregate(depth~chrom+pos+ref+alt, pd, min)
        min.depth <- min.depth[min.depth$depth >= input$depth,
                               c("chrom", "pos", "ref", "alt")]
        pd <- merge(pd, min.depth, by=c("chrom", "pos", "ref", "alt"),
                    suffixes=c("", ".min"))

        min.freqs <- aggregate(vaf~chrom+pos+ref+alt, pd, min)
        min.freqs <- subset(min.freqs, vaf <= input$maxmin)
        pd <- merge(pd, min.freqs, by=c("chrom", "pos", "ref", "alt"),
                    suffixes=c("", ".min"))

        if (!input$all.chrom) {
            pd <- subset(pd, chrom %in% input$chrom)
        }
        
        if (input$hide.silent) {
            pd <- subset(pd, class != "Silent")
        }

        if (nrow(pd) == 0) return (add.zero.row(pd))

        max.date <- max(pd$date)
        late.samples <- unique(subset(pd, is.na(date) | date == max.date)$sample)
        pd <- subset(pd, sample == input$sample | !sample %in% late.samples)

        if (nrow(pd) == 0) return(add.zero.row(pd))

	    if (input$order == "Highest fraction") {
	        aggfun <- sum
	    } else if (input$order == "Most change") {
	        aggfun <- function (x) max(abs(diff(x)))
	    }
	    agg <- aggregate(vaf~chrom+pos+ref+alt, pd, aggfun)
	    agg <- agg[order(-agg$vaf),]
	    agg <- head(agg, input$n)
	    pd <- merge(pd, agg, by=c("chrom", "pos", "ref", "alt"), 
              suffixes=c("", ".agg"))

        if (nrow(pd) == 0) return (add.zero.row(pd))
        droplevels(pd[order(pd$date),])
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

    output$chrom.select <- renderUI({
        html <- selectInput("chrom", "Chromosome:", levels(d$chrom), multi=T, selected=input$chrom)
        if (input$all.chrom) 
            gsub("(<select[^>]*)>", "\\1 disabled>", html)
        else
            html
    })

    output$sample.select <- renderUI({
        s <- subset(d, patient==input$patient, select=c("sample", "date"))
        max.date <- max(s$date)
        choices <- unique(subset(s, is.na(date) | date == max.date)$sample)
        html <- selectInput("sample", "Late sample:", as.character(choices))
        if (length(choices) == 1)
            gsub("(<select[^>]*)>", "\\1 disabled>", html)
        else
            html
    })

    freqPlot.vis %>% bind_shiny("freqPlot")

    output$freqTable <- renderDataTable({
        agg
    })
})
