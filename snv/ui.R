library(shiny)
library(ggvis)

shinyUI(navbarPage("VAF",

  tabPanel("By patient",
    sidebarLayout(
      sidebarPanel(
        uiOutput("patient.select"),
        uiOutput("sample.select"),
        uiOutput("chrom.select"),
        checkboxInput("all.chrom", "Select all chromosomes", value=T),
        br(),
        sliderInput("depth", "Minimum coverage depth:", min=0, max=100, step=1, value=0),
        selectInput("color", "Color by:", c("None", "Chromosome", "Mutation type")),
        checkboxInput("hide.silent", "Hide silent mutations"),
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

  tabPanel("Overall",
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
        fluidRow( column(12, dataTableOutput("freqTable") ))
    )
  )
))
