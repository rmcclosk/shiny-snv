library(shiny)
library(ggplot2)
library(ggdendro)

shinyServer(function(input, output) {
    plot.samples <- reactive({
        m <- subset(metadata, patient==input$patient & time.point > 0)
        sapply(sort(unique(m$time.point)), function (tp) {
            s <- input[[paste0("tp", tp)]]
            if (is.null(s))
                subset(m, time.point==tp)$sample[1]
            else
                s
        })
    })

    plot.purity <- reactive({
        m <- subset(metadata, sample %in% plot.samples())
        d <- data.frame(x=c(rep(0, nrow(m)), rep(max(plot.segments()$pos), nrow(m))), 
                        y=rep(m$purity/100, 2),
                        sample=rep(as.character(m$sample), 2))
        d[order(d$sample),]
    })

    plot.segments <- reactive({
        ps <- subset(segments, patient==input$patient)
        ps$sample <- as.character(ps$sample)
        ps <- data.frame(pos=c(ps$adj.start, ps$adj.end), 
                         copy.number=c(ps$copy.number, ps$copy.number),
                         sample=c(ps$sample, ps$sample),
                         segment=rep(1:nrow(ps), 2),
                         time.point=c(ps$time.point, ps$time.point))

        samples <- plot.samples()

        adj <- (1:length(samples))-floor(length(samples)/2)
        tmp <- data.frame(sample=samples, adj=adj)

        ps <- ps[ps$sample %in% samples,]
        ps <- merge(ps, tmp)
        ps$copy.number <- ps$copy.number + 0.1*ps$adj
        droplevels(ps[order(ps$sample, ps$segment),])
    })

    plot.variants <- reactive({
        pv <- subset(variants, patient==input$patient & chrom==input$chrom)
        pv$sample <- as.character(pv$sample)
        pv <- data.frame(pos=c(pv$pos, pv$pos),
                         vaf=c(pv$vaf, pv$vaf.corrected),
                         sample=c(pv$sample, pv$sample),
                         purity=c(pv$purity/100, pv$purity/100))

        samples <- plot.samples()

        adj <- (1:length(samples))-floor(length(samples)/2)
        tmp <- data.frame(sample=samples, adj=adj)

        pv <- pv[pv$sample %in% samples,]
        pv <- merge(pv, tmp)
        if (nrow(pv) > 0)
            pv$pos <- pv$pos + max(pv$pos)/200*pv$adj
        droplevels(pv)
    })

    output$sample.select <- renderUI({
        m <- subset(metadata, patient==input$patient & time.point != 0)
        do.call(tagList, as.list(by(m, m$time.point, function (s) {
            tp <- s[1, "time.point"]
            var <- paste0("tp", tp)
            label <- paste0("Time point ", tp, " sample:") 
            choices <- unique(as.character(s$sample))
            selectInput(var, label, choices)
        })))
    })

    tooltip <- function (x) {
        pd <- isolate(plot.variants())
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

    segPlot.vis <- reactive({
        plot.segments %>% 
            ggvis(x=~pos, y=~copy.number) %>%
            add_axis("x", title="position") %>%
            add_axis("y", title="copy number") %>%
            scale_numeric("y", domain=c(0, 5), nice=F, round=T) %>%
            group_by(segment) %>%
            layer_paths(stroke=~sample, fill=~sample, strokeWidth:=10)
    })

    segPlot.vis %>% bind_shiny("segPlot")

    output$hclust <- renderPlot({
        n.samples <- nrow(subset(metadata, patient==input$patient & time.point > 0))
        if (n.samples > 2)
            ggdendrogram(hclusts[[input$patient]])
    })
})
