library(shiny)
library(ggplot2)
library(ggdendro)

shinyServer(function(input, output) {

    # get samples which will be displayed
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

    # get segments to plot
    plot.segments <- reactive({
        ps <- subset(segments, patient==input$patient)
        ps$sample <- as.character(ps$sample)
        samples <- plot.samples()

        adj <- (1:length(samples))-floor(length(samples)/2)
        tmp <- data.frame(sample=samples, adj=adj)

        ps <- ps[ps$sample %in% samples,]
        ps <- merge(ps, tmp)
        ps$copy.number <- ps$copy.number + 0.1*ps$adj
        ps <- merge(ps, aggregate(prevalence~sample, ps, max))
        droplevels(ps)
    })

    # get SNV's to plot
    plot.variants <- reactive({
        d <- subset(variants, sample %in% plot.samples())
        d$vaf.corrected[d$vaf.corrected > 1.5] <- 1.5
        d
    })

    # render select boxes for patient samples
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

    # plot segments
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

    output$vafPlot <- renderPlot({
        d <- plot.variants()
        starts <- sort(unique(d$chr.start))
        breaks <- starts + c(diff(starts)/2, 50000000)
        ggplot(d, aes(x=adj.pos, y=vaf.corrected, color=sample)) +
            geom_point() +
            theme_bw() +
            ylab("corrected VAF") +
            xlab("chromosome") +
            theme(axis.ticks.x=element_blank()) +
            geom_segment(aes(x=chr.start, xend=chr.start, y=0, yend=1.5, group=chrom), color="grey", linetype="dashed") +
            scale_x_continuous(breaks=breaks, labels=c(1:22, "X"))
    })

    # plot hierarchical clustering
    output$hclust <- renderPlot({
        n.samples <- nrow(subset(metadata, patient==input$patient & time.point > 0))
        if (n.samples > 2)
            ggdendrogram(hclusts[[input$patient]])
    })
})
