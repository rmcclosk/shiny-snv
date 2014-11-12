library(shiny)
library(ggplot2)

shinyServer(function(input, output) {
    plot.segments <- reactive({
        ps <- subset(segments, patient==input$patient & chrom==input$chrom)
        ps$sample <- as.character(ps$sample)
        ps <- data.frame(pos=c(ps$start, ps$end), 
                         copy.number=c(ps$copy.number, ps$copy.number),
                         sample=c(ps$sample, ps$sample),
                         segment=rep(1:nrow(ps), 2),
                         time.point=c(ps$time.point, ps$time.point))

        samples <- sapply(unique(ps$time.point), function (tp) {
            input[[paste0("tp", tp)]]
        })

        ps <- ps[ps$sample %in% samples,]
        droplevels(ps[order(ps$sample, ps$segment),])
    })

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

    segPlot.vis <- reactive({
        plot.segments %>% 
            ggvis(x=~pos, y=~copy.number) %>%
            add_axis("x", title="position") %>%
            add_axis("y", title="copy number") %>%
            scale_numeric("y", domain=c(0, 4), nice=F, round=T) %>%
            group_by(segment) %>%
            layer_paths(stroke=~sample, fill=~sample, strokeWidth:=10)
    })

    segPlot.vis %>% bind_shiny("segPlot")
})
