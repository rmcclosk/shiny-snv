library(shiny)
library(ggvis)

shinyUI(navbarPage("",

  tabPanel("VAF by patient",
    sidebarLayout(
      sidebarPanel(
        selectInput("patient.vaf", "Patient:", unique(variants$patient)),
        uiOutput("sample.select.vaf"),
        uiOutput("chr.select"),
        checkboxInput("all.chr", "Select all chromosomes", value=T),
        uiOutput("gene.select"),
        checkboxInput("all.gene", "Select all genes", value=T),
        br(),
        sliderInput("depth", "Minimum coverage depth:", min=0, max=100, step=1, value=0),
        selectInput("color", "Color by:", c("None", "Chromosome", "Mutation type")),
        checkboxInput("hide.silent.plot", "Hide silent mutations"),
        br(),
        sliderInput("n", "Number to show:", min=0, max=100, step=1, value=20),
        selectInput("order", "Order by:", c("Highest fraction", "Most change")),
        sliderInput("maxmin", "Must fall below:", min=0, max=1, value=1),
        checkboxInput("correct.purity", "Correct VAF to CCF")
      ),
  
      mainPanel(
        ggvisOutput("freqPlot")
      )
    )
  ),

  tabPanel("VAF overall",
    fluidPage(
        fluidRow( 
            wellPanel(
                column(10,
                    checkboxInput("hide.silent.table", "Hide silent mutations")
                ),
                column(2,
                    downloadButton('downloadData', 'Download')
                )
            )
        ),
        fluidRow( column(12, dataTableOutput("table") ))
    )
  ),

  tabPanel("CNA by patient",
    sidebarLayout(
      sidebarPanel(
        selectInput("patient.cna", "Patient:", unique(segments$patient)),
        uiOutput("sample.select.cna")
      ),
  
      mainPanel(
        plotOutput("segsPlot")
      )
    )
  )
))
