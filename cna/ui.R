library(shiny)

shinyUI(fluidPage(
  titlePanel("Copy number alterations"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient", "Patient ID:", levels(d$patient), multi=F),
      uiOutput("sample.select"),
      selectInput("chrom", "Chromosome:", c(1:22, "X"))
    ),

    mainPanel(
      uiOutput("plots")
    )
  )
)
