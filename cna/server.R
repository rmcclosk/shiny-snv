library(shiny)
library(ggplot2)

shinyServer(function(input, output) {
    output$segments <- renderImage({
        list(src = file.path("plots", paste0(input$patient, ".png")),
             contentType = "image/png")
    }, deleteFile=F)
})
