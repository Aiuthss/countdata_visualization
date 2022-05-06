#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

shinyUI(
  fluidPage(
    sidebarPanel(width=3, 
      selectInput(inputId="strain", label="select strain", choices = c("MMusculus", "HSapiens", "MAuratus")),
      fileInput(inputId="fileInput", label="Upload file", multiple=FALSE),
      fileInput(inputId="metadata_file", label="Upload metadata file", multiple=FALSE),
      radioButtons(inputId="grouping", label="grouping", choices=c("manual", "metadata file")),
      actionButton(inputId = "confirm_button", label="confirm"),
    ),
    mainPanel(
      tabsetPanel(type="tabs",
        tabPanel("table",
          textOutput(outputId = "textOutput"),
          dataTableOutput(outputId = "table"),
          dataTableOutput(outputId = "test_table")
        ),
        tabPanel("analysis",
          dataTableOutput(outputId = "res")
        ),
        tabPanel("MAplot",
          fluidRow(
            column(8,
              wellPanel(
                textOutput(outputId = "plot_title"),
                plotOutput(outputId = "plot1"),
                plotOutput(outputId = "plot2"),
                plotOutput(outputId = "vol_plot1"),
                plotOutput(outputId = "vol_plot2")
              )
            ),
            column(4,
              wellPanel(
                uiOutput("grouping_input"),
                uiOutput("plot_parameter")
              )
            )
          )
        ),
        tabPanel("GSEA",
          selectInput(inputId="ontology", label="select ontology", choices = c("Biological Process", "Cellular Component", "Molecular Function")),
          actionButton(inputId="GSEA_button", label="run"),
          tabsetPanel(type="tabs",
            tabPanel("Upregulated",
              dataTableOutput(outputId = "go_up")
            ),
            tabPanel("Downregulated",
              dataTableOutput(outputId = "go_down")
            )
          )
        )
      )
    )
  )
)