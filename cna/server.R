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

    plot.segments <- reactive({
        ps <- subset(segments, patient==input$patient)
        ps$sample <- as.character(ps$sample)
        #ps <- data.frame(pos=c(ps$adj.start, ps$adj.end), 
        #                 copy.number=c(ps$copy.number, ps$copy.number),
        #                 sample=c(ps$sample, ps$sample),
        #                 segment=rep(1:nrow(ps), 2),
        #                 time.point=c(ps$time.point, ps$time.point))

        samples <- plot.samples()

        adj <- (1:length(samples))-floor(length(samples)/2)
        tmp <- data.frame(sample=samples, adj=adj)

        ps <- ps[ps$sample %in% samples,]
        ps <- merge(ps, tmp)
        ps$copy.number <- ps$copy.number + 0.1*ps$adj
        #droplevels(ps[order(ps$sample, ps$segment),])
        droplevels(ps)
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

    output$segPlot <- renderPlot({
        d <- plot.segments()
        starts <- sort(unique(d$chr.start))
        breaks <- starts + c(diff(starts)/2, 50000000)
        ggplot(d, aes(x=adj.start, y=copy.number, color=sample)) +
            geom_segment(aes(xend=adj.end, yend=copy.number), size=3) +
            theme_bw() +
            ylab("copy number") +
            xlab("chromosome") +
            theme(axis.ticks.x=element_blank()) +
            geom_segment(aes(x=chr.start, xend=chr.start, y=0, yend=5, group=chrom), color="grey", linetype="dashed") +
            scale_x_continuous(breaks=breaks, labels=c(1:22, "X"))
    })

#    segPlot.vis <- reactive({
#        plot.segments %>% 
#            ggvis(x=~pos, y=~copy.number) %>%
#            add_axis("x", title="position") %>%
#            add_axis("y", title="copy number") %>%
#            scale_numeric("y", domain=c(0, 5), nice=F, round=T) %>%
#            group_by(segment) %>%
#            layer_paths(stroke=~sample, fill=~sample, strokeWidth:=10) %>%
#            layer_paths(data=chr.bounds, x=~x, y=~y, strokeWidth=~y) %>%
#            scale_numeric("strokeWidth", range=c(0, 2))
#    })
#
#    segPlot.vis %>% bind_shiny("segPlot")

    output$hclust <- renderPlot({
        n.samples <- nrow(subset(metadata, patient==input$patient & time.point > 0))
        if (n.samples > 2)
            ggdendrogram(hclusts[[input$patient]])
    })
})
