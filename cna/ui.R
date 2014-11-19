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
      #ggvisOutput("segPlot"),
      plotOutput("segPlot"),
      br(),
      plotOutput("hclust")
    )
  )
))
