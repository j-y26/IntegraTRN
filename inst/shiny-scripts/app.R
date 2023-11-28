# Purpose: Define the shiny app for the package
# Author: Jielin Yang
# Date: 2023-11-29
# Version: 1.0
# Bugs and Issues: None

# Load libraries
library(shiny)
library(shinyalert)

# Define the UI for the application
ui <- fluidPage(

  # Title of the application, consistent with Github repo
  titlePanel(tags$h1(tags$b("IntegraTRN:"), "Integrating multi-omics data for 
    the exploration of regulatory mechanisms and the inference of core TRNs 
    underlying transcriptomic alterations")),
  
  # Use sidebar layout with a sidebar panel and a main panel
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Description of the app
      tags$p(tags$b("Description:"), "The IntegraTRN shiny app performs 
        integrative analysis of multi-omics data to infer the core
        transcriptional regulatory networks (TRNs) that explains the
        condition-specific transcriptomic alterations. The app takes
        as inputs a set of expression data (RNAseq, small RNAseq, and/or
        proteomics), a set of chromatin accessible regions (ATACseq), and
        a set of externally curated regulatory interactions, such as TF-target
        gene and miRNA-target gene interactions. The app then performs
        integrative analysis to infer the core TRNs, including predictive
        modeling with a random forest or extra-trees classifier, to predict
        the regulatory interactions, with rigorous filtering to preserve
        high-confidence interactions."),

      # Empty space
      br(),
      
      # Instructions for the app
      tags$b("Below, please upload the input data files tailored
              to each data type according to the detailed instructions
              below. Then press 'Run'. Navigate through the different tabs
              to the right to explore the results. According to the data type
              provided, the app will automatically perform the corresponding
              analysis. Not all tabs will be available if the corresponding
              data type is not provided."),
      
      # Empty space
      br(),
      br(),
      br(),

      # Data inputs

      # RNAseq data input

      # RNAseq count matrix
      tags$b("Upload RNAseq count matrix (.csv):"),
      fileInput(inputId = "rnaseq_count_matrix",
                label = "A numeric matrix with genes as rows 
                and samples as columns. The gene names should be specified in 
                the first column. The first row should be sample names.",
                accept = c(".csv")),
      div(style = "margin-top: -35px"),    # move the download link up
      # The above solution to adjust the position of the download link is
      # from https://stackoverflow.com/questions/62476054/reduce-space-between-
      # fileinput-and-text-in-shiny by Stéphane Laurent
      downloadLink(outputId = "download_rnaseq_count_matrix",
                   label = "Download example RNAseq count matrix for human
                            fetal heart development (.csv)"),
      br(),
      br(),

      # RNAseq sample metadata
      tags$b("Upload RNAseq sample metadata (.csv):"),
      fileInput(inputId = "rnaseq_sample_metadata",
                label = "A data frame with samples as rows and attributes
                as columns. The first column should be sample names. The
                first row should be metadata names.",
                accept = c(".csv")),
      div(style = "margin-top: -35px"),   # move the download link up
      downloadLink(outputId = "download_rnaseq_sample_metadata",
                   label = "Download example RNAseq sample metadata for human
                            fetal heart development (.csv)"),
      br(),
      # Select the grouping variable for RNAseq samples
      uiOutput("rnaseq_sample_attributes")
    ),

    # Main panel for displaying outputs
    mainPanel(

    ),

  ),
)


# Define the server for the application
server <- function(input, output) {
  
  # RNAseq data input

  # RNAseq count matrix
  rnaseq_count_matrix <- reactive({
    req(input$rnaseq_count_matrix)

    # Enforce correct file extension
    ext <- tools::file_ext(input$rnaseq_count_matrix$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))

    rnaseq <- read.csv(input$rnaseq_count_matrix$datapath, 
                       header = TRUE, 
                       row.names = 1)
    rnaseq <- as.matrix(rnaseq)
    return(rnaseq)
  })
  # for downloading example RNAseq count matrix
  output$download_rnaseq_count_matrix <- downloadHandler(
    filename = function() {
      paste("rnaseq_count_matrix_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(RNAseq_heart, file, row.names = TRUE)
    }
  )

  # RNAseq sample metadata
  rnaseq_sample_metadata <- reactive({
    req(input$rnaseq_sample_metadata)

    # Enforce correct file extension
    ext <- tools::file_ext(input$rnaseq_sample_metadata$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))

    read.csv(input$rnaseq_sample_metadata$datapath, 
             header = TRUE, 
             row.names = 1)
  })
  # for downloading example RNAseq sample metadata
  output$download_rnaseq_sample_metadata <- downloadHandler(
    filename = function() {
      paste("rnaseq_sample_metadata_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(RNAseq_heart_samples, file, row.names = TRUE)
    }
  )
  # Retrieve the attribute (column) names of the sample metadata
  rnaseq_sample_attributes <- reactive({
    req(input$rnaseq_sample_metadata)
    colnames(rnaseq_sample_metadata())
  })
  # Render selector UI for sample attributes for RNAseq
  output$rnaseq_sample_attributes <- renderUI({
    selectInput(inputId = "rnaseq_sample_grouping",
                label = "Select the grouping variable for RNAseq samples:",
                choices = rnaseq_sample_attributes(),
                multiple = FALSE)
  })



  
}



# Build the shiny app object
shinyApp(ui = ui, server = server)