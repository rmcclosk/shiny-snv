library(shiny)
library(ggvis)

shinyUI(fluidPage(
  titlePanel("Variant allele fractions through time"),

  sidebarLayout(
    sidebarPanel(
      selectInput("patient", "Patient ID:", levels(d$patient), multi=F),
      uiOutput("sample.select"),
      uiOutput("chrom.select"),
      checkboxInput("all.chrom", "Select all chromosomes"),
      br(),
      sliderInput("depth", "Minimum coverage depth:", min=0, max=100, step=1, value=0),
      selectInput("color", "Color by:", c("None", "Chromosome", "Mutation type")),
      checkboxInput("hide.silent", "Hide silent mutations"),
      br(),
      sliderInput("n", "Number to show:", min=0, max=100, step=1, value=20),
      selectInput("order", "Order by:", c("Highest fraction", "Most change")),
      sliderInput("maxmin", "Must fall below:", min=0, max=1, value=1)
    ),

    mainPanel(
      ggvisOutput("freqPlot")
    )
  )
))
