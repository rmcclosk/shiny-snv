library(shiny)
library(ggvis)

shinyUI(fluidPage(
  titlePanel("Copy number alterations"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient", "Patient ID:", levels(segments$patient)),
      uiOutput("sample.select")
    ),

    mainPanel(
      plotOutput("segPlot"),
      br(),
      plotOutput("vafPlot"),
      br(),
      plotOutput("hclust")
    )
  )
))
