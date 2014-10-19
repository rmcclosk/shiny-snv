library(shiny)

shinyUI(fluidPage(
  titlePanel("Copy number alterations"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient", "Patient ID:", levels(sample.data$patient))
    ),

    mainPanel(
      imageOutput("segments")
    )
  )
))
