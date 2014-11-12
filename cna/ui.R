library(shiny)
library(ggvis)

shinyUI(fluidPage(
  titlePanel("Copy number alterations"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient", "Patient ID:", levels(segments$patient)),
      uiOutput("sample.select"),
      selectInput("chrom", "Chromosome:", levels(segments$chrom))
    ),

    mainPanel(
      ggvisOutput("segPlot")
    )
  )
))
