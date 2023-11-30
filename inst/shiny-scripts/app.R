# Purpose: Define the shiny app for the package
# Author: Jielin Yang
# Date: 2023-11-29
# Version: 1.0
# Bugs and Issues: None

# Load libraries
library(shiny)
library(shinyalert)

# Define some global variables
# Note that RNA, SMALLRNA, ATAC, and PROTEIN is already defined within the
# package
# The following variables are used to define external data types
# Note these definitions are only available within the app
MIRNA <- "miRNA_target"
TF <- "TF_target"


# Define the UI for the application
ui <- fluidPage(
  # Since fileInput labels are detailed instructions, custom CSS
  # is used to make the labels not bold
  tags$style(HTML(".control-label { font-weight: normal; }")),

  # Also since there are multi-level headers, to better format the text,
  # custom CSS is used to define bold and italic headers
  tags$style(HTML("h5 { font-weight: bold; 
                        font-style: italic; 
                        font-size: 1.1em; }")),
  
  tags$h4("Section 1: Raw Data Input"),

  # Title of the application, consistent with Github repo
  titlePanel(tags$h1(tags$b("IntegraTRN:"), "Integrating multi-omics data for 
    the exploration of regulatory mechanisms and the inference of core TRNs 
    underlying transcriptomic alterations")),
  
  # Use sidebar layout with a sidebar panel and a main panel
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(

      # Creating multiple tabs for different input sections
      tabsetPanel(

        # === Panel for introduction ===========================================
        tabPanel("Introduction",
        # Panel title
        tags$h4("Introduction"),
        br(),

        # Description of the app
        tags$b("Package Description:"),
        tags$p("The IntegraTRN shiny app 
          performs integrative analysis of multi-omics data to infer 
          the core transcriptional regulatory networks (TRNs) that 
          explains the condition-specific transcriptomic alterations. 
          The app takes as inputs a set of expression data (RNAseq, 
          small RNAseq, and/or proteomics), a set of chromatin 
          accessible regions (ATACseq), and a set of externally 
          curated regulatory interactions, such as TF-target
          gene and miRNA-target gene interactions. The app then 
          performs integrative analysis to infer the core TRNs, 
          including predictive modeling with a random forest or 
          extra-trees classifier, to predict the regulatory 
          interactions, with rigorous filtering to preserve
          high-confidence interactions."),

        # Empty space
        br(),
        
        # Image summary of the package
        tags$img(src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/Schematics.jpg",   # link cannot be broken to shorten the line
                  width = "550px", align = "center"),
        tags$p("Schematics of the IntegraTRN package workflow"),

        br(),
        br(),

        # Instructions for the app
        tags$b("Instructions:"),
        tags$b("In the following sections, please upload the input 
                data files tailored to each data type according to the 
                detailed instructions. Then press 'Run'. Navigate 
                through the different tabs to the right to explore 
                the results. According to the data type provided, 
                the app will automatically perform the corresponding
                analysis. Not all tabs will be available if the 
                corresponding data type is not provided."),
        ),

        # === Panel for section 1: select omics data ===========================
        tabPanel("Section 1",
        tags$h4("Section 1: Select Omics Data"),
        br(),

        # Description of the section
        tags$p(tags$b("Description:"), "In this section, please select the 
          omics data types that you would like to analyze. The app will 
          automatically perform the corresponding analysis based on the 
          data types provided."),
        br(),
        tags$p(tags$b("Note:"), "RNAseq data is required for the analysis. 
          Additional data types are optional. However, it is recommended
          to provide a comprehensive set of data to generate a more biologically
          relevant TRN. Example data sets for each type are provided for 
          download in the following sections."),
        br(),
        tags$p(tags$b("Note:"), "Select the type of data for the analysis by 
          checking the corresponding boxes. Hover over the text or check the 
          boxes to learn about each data type."),

        # For each type, provide an image on the left and a checkbox with the
        # name of the type on the right
        # RNAseq (not a checkbox, required)
        column(width = 6,
               tags$img(src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_rnaseq.jpg",
                        width = "200px", align = "center"),
               div(style = "margin-top: 10px"),
               tags$b(id = "omicRNA", "RNAseq (required)"), align = "center"),
               bsTooltip(id = "omicRNA",
                         title = paste0("Bulk RNA sequencing (RNAseq) data is ",
                         "required for the analysis, which provides a full ",
                         "transcriptomic profile of gene expression for the ",
                         "tissue/biological system of interest.")),
        # The use of bsTooltip is to add hover-over tooltips is derived from
        # https://stackoverflow.com/questions/16449252/tooltip-on-shiny-r
        # by Praneeth Jagarapu

        # Small RNAseq
        column(width = 6,
               tags$img(src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_smallrnaseq.jpg",
                        width = "200px", align = "center"),
               checkboxInput(inputId = "omicSmallRNA", 
                             label = tags$b("Small RNAseq"),
                             value = FALSE), align = "center"),
               bsTooltip(id = "omicSmallRNA",
                         title = paste0("Small RNA sequencing (small RNAseq) ",
                          "data is optional for the analysis, which provides ",
                          "insights into the expression of small RNAs, such ",
                          "as miRNAs, piRNAs, and snoRNAs. They are important ",
                          "regulators of gene expression.")),
        
        # Proteomics
        column(width = 6,
               tags$img(src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_proteomics.jpg",
                        width = "200px", align = "center"),
               checkboxInput(inputId = "omicProteomics", 
                             label = tags$b("Proteomics"),
                             value = FALSE), align = "center"),
               bsTooltip(id = "omicProteomics",
                          title = paste0("Proteomics data is optional for the ",
                            "analysis. Cross referencing transcriptomic and ",
                            "proteomic data can further validate gene ",
                            "expression changes at the protein level.")),
        
        # ATACseq
        column(width = 6,
               tags$img(src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_atac.jpg",
                        width = "200px", align = "center"),
               checkboxInput(inputId = "omicATAC", 
                             label = tags$b("ATACseq"),
                             value = FALSE), align = "center"),
               bsTooltip(id = "omicATAC",
                          title = paste0("ATACseq data is optional. It ",
                            "provides valuable information on chromatin ",
                            "accessibility, which reflects putative ",
                            "binding sites of transcription factors. ")),
        
        # External data
        # miRNA - target gene interactions
        column(width = 6,
               tags$img(src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_mirna.jpg",
                        width = "200px", align = "center"),
               checkboxInput(inputId = "externalmiRNA", 
                             label = tags$b("microRNA - Target Interactions"),
                             value = FALSE), align = "center"),
               bsTooltip(id = "externalmiRNA",
                          title = paste0("microRNAs (miRNAs) are important ",
                            "regulators of gene expression by binding to ",
                            "target genes and leading to their degradation ",
                            "or translational repression. Curated miRNA - ",
                            "target gene interactions are available from ",
                            "external databases, and can be used to establish ",
                            "a non-tissue-specific context for the analysis.")),
        
        # TF - target gene interactions
        column(width = 6,
               tags$img(src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_tf.jpg",
                        width = "200px", align = "center"),
               checkboxInput(inputId = "externalTF", 
                             label = tags$b("TF - Target Interactions"),
                             value = FALSE), align = "center"),
               bsTooltip(id = "externalTF",
                          title = paste0("Transcription factors (TFs) ",
                            "can mediate transcriptional activation or ",
                            "repression of target genes. Curated TF - ",
                            "target gene interactions are available from ",
                            "previous ChIP-seq studies, but may lack ",
                            "tissue-specificity. Intersecting with ATACseq ",
                            "data can help to identify tissue-specific ",
                            "TF - target gene interactions.")),

        br(),
        tags$p("Please upload the selected data types in the following 
                sections."),
        br(),
        
        tags$h5("Example Data Sets"),

        tags$p("If you would like to take a look at the example data sets, the
                following links provide a centralized access to the data sets.
                Example data are also available for download in the following
                sections."),
        tags$p("Example datasets for human fetal heart development:"),

        tags$ul(
          tags$li(downloadLink(outputId = "downloadRnaseqCountMatrix", 
                               label = "RNAseq count matrix (.csv)")),
          tags$li(downloadLink(outputId = "downloadRnaseqSampleMetadata", 
                               label = "RNAseq sample metadata (.csv)")),
          tags$li(downloadLink(outputId = "downloadSmallRnaseqCountMatrix", 
                               label = "Small RNAseq count matrix (.csv)")),
          tags$li(downloadLink(outputId = "downloadSmallRnaseqSampleMetadata", 
                               label = "Small RNAseq sample metadata (.csv)")),
          tags$li(downloadLink(outputId = "downloadProteomicsCountMatrix", 
                               label = "Proteomics count matrix (.csv)")),
          tags$li(downloadLink(outputId = "downloadProteomicsSampleMetadata", 
                               label = "Proteomics sample metadata (.csv)")),
          tags$li(downloadLink(outputId = "downloadATACPeak1", 
                               label = "ATACseq peak for condition 1 (.bed)")),
          tags$li(downloadLink(outputId = "downloadATACPeak2", 
                               label = "ATACseq peak for condition 2 (.bed)")),
          tags$li(downloadLink(outputId = "downloadmiRNATarget", 
                               label = "miRNA - target interactions (.csv)")),
          tags$li(downloadLink(outputId = "downloadTFTarget", 
                               label = "TF - target interactions (.csv)")),
          tags$li(downloadLink(outputId = "downloadsmallRNAAnnotation", 
                               label = "Small RNA annotation (.csv)")),
          tags$li(downloadLink(outputId = "downloadproteinNameConversion", 
                               label = "Gene - protein name conversion (.csv)"))
        ),
      ),



      #   # === Panel for section 2: raw data input ==============================
      #   tabPanel("Section 2",
      #   tags$h4("Section 2: Raw Data Input"),
      #   br(),

      #   # Description of the section
      #   tags$p(tags$b("Description:"), "In this section, please upload 
      #     the raw data files for each data type."),
      #   br(),

      #   # Data inputs

      #   # RNAseq data input
      #   tags$h5("Part 1: RNAseq Data Input (Required)"),

      #   # RNAseq count matrix
      #   tags$b("Upload RNAseq count matrix (.csv):"),
      #   fileInput(inputId = "rnaseqCountMatrix",
      #             label = "A numeric matrix with genes as rows 
      #             and samples as columns. The gene names should be specified in 
      #             the first column. The first row should be sample names.",
      #             accept = c(".csv")),
      #   div(style = "margin-top: -35px"),    # move the download link up
      #   # The above solution to adjust the position of the download link is
      #   # from https://stackoverflow.com/questions/62476054/reduce-space-between-
      #   # fileinput-and-text-in-shiny by Stéphane Laurent
      #   downloadLink(outputId = "downloadRnaseqCountMatrix",
      #               label = "Download example RNAseq count matrix for human
      #                         fetal heart development (.csv)"),
      #   br(),
      #   br(),

      #   # RNAseq sample metadata
      #   tags$b("Upload RNAseq sample metadata (.csv):"),
      #   fileInput(inputId = "rnaseqSampleMetadata",
      #             label = "A data frame with samples as rows and attributes
      #             as columns. The first column should be sample names. The
      #             first row should be metadata names.",
      #             accept = c(".csv")),
      #   div(style = "margin-top: -35px"),   # move the download link up
      #   downloadLink(outputId = "downloadRnaseqSampleMetadata",
      #               label = "Download example RNAseq sample metadata for human
      #                         fetal heart development (.csv)"),
      #   br(),
      #   # Select the grouping variable for RNAseq samples
      #   uiOutput("rnaseqSampleGrouping"),

      #   br(),
      #   br(),

      #   # Small RNAseq data input
      #   tags$h5("Part 2: Small RNAseq Data Input (Optional)"),

      #   # Small RNAseq count matrix
      #   tags$b("Upload small RNAseq count matrix (.csv):"),
      #   fileInput(inputId = "smallRnaseqCountMatrix",
      #             label = "A numeric matrix with transcripts as rows 
      #             and samples as columns. The small RNA names should be specified
      #             in the first column. The first row should be sample names.",
      #             accept = c(".csv")),
      #   div(style = "margin-top: -35px"),    # move the download link up
      #   downloadLink(outputId = "downloadSmallRnaseqCountMatrix",
      #                label = "Download example small RNAseq count matrix for
      #                         human fetal heart development (.csv)"),
      #   br(),
      #   br(),

      #   # Small RNAseq sample metadata
      #   tags$b("Upload small RNAseq sample metadata (.csv):"),
      #   fileInput(inputId = "smallRnaseqSampleMetadata",
      #             label = "A data frame with samples as rows and attributes
      #             as columns. The first column should be sample names. The
      #             first row should be metadata names.",
      #             accept = c(".csv")),
      #   div(style = "margin-top: -35px"),   # move the download link up
      #   downloadLink(outputId = "downloadSmallRnaseqSampleMetadata",
      #               label = "Download example small RNAseq sample metadata for
      #                         human fetal heart development (.csv)"),
      #   br(),
      #   # Select the grouping variable for small RNAseq samples
      #   uiOutput("smallRnaSampleGrouping"),

      #   br(),
      #   br(),

      #   # Proteomics data input
      #   tags$h5("Part 3: Proteomics Data Input (Optional)"),

      #   # Proteomics count matrix
      #   tags$b("Upload proteomics count matrix (.csv):"),
      #   fileInput(inputId = "proteomicsCountMatrix",
      #             label = "A numeric matrix with proteins as rows 
      #             and samples as columns. The protein names should be specified
      #             in the first column. The first row should be sample names.",
      #             accept = c(".csv")),
      #   div(style = "margin-top: -35px"),    # move the download link up
      #   downloadLink(outputId = "downloadProteomicsCountMatrix",
      #               label = "Download example proteomics count matrix for
      #                         human fetal heart development (.csv)"),
      #   br(),
      #   br(),

      #   # Proteomics sample metadata
      #   tags$b("Upload proteomics sample metadata (.csv):"),
      #   fileInput(inputId = "proteomicsSampleMetadata",
      #             label = "A data frame with samples as rows and attributes
      #             as columns. The first column should be sample names. The
      #             first row should be metadata names.",
      #             accept = c(".csv")),
      #   div(style = "margin-top: -35px"),   # move the download link up
      #   downloadLink(outputId = "downloadProteomicsSampleMetadata",
      #               label = "Download example proteomics sample metadata for
      #                         human fetal heart development (.csv)"),
      #   br(),
      #   # Select the grouping variable for proteomics samples
      #   uiOutput("proteomicsSampleGrouping"),
      # ),

      #   # === Panel for section 3: annotations and external data ===============
      #   tabPanel("Section 3",
      #   tags$h4("Section 3: Annotations and Externally Curated Regulatory
      #           Interactions"),
      #   br(),

      #   # Description of the section
      #   tags$p(tags$b("Description:"), "In this section, please upload 
      #     annotation files for small RNAs and proteins, as well as
      #     externally curated regulatory interactions, such as TF-target
      #     gene and miRNA-target gene interactions."),
      #   br(),

      #   tags$p("An example tool to search for externally curated interactions is
      #           miRNet, available ",
      #           tags$a(href = "https://www.mirnet.ca/", "miRNet"), "."),

        

      #   # Space holder, implement later

      #   ),

      #   # === Panel for section 4: integrative analysis ========================
      #   tabPanel("Section 4",
      #   tags$h4("Section 4: Integrative Analysis"),
      #   br(),

      #   # Description of the section
      #   tags$p(tags$b("Description:"), "In this section, the app performs
      #     integrative analysis of the multi-omics data to infer the core
      #     TRNs that explains the condition-specific transcriptomic
      #     alterations. Please provide settings to the program to define
      #     the behavior of the analysis."),
      #   br(),

      #   # Space holder, implement later

      #   ),
      ),
    ),

    # Main panel for displaying outputs
    mainPanel(
      # Create a tabset panel for different output sections
      tabsetPanel(
        




        # ======================================================================
        # === FOR OUTPUT TESTING ONLY, REMOVE LATER ============================
        # ======================================================================
        tabPanel("Output Testing",
        tags$h4("Output Testing"),
        br(),

        # Test correct selection of omics data types
        tags$p("The following text should be the selected omics data types:"),
        textOutput(outputId = "omicsDataTypes"),
        br(),










        ),
      
      
      
      
      
      ),
    ),
  ),
  
  # Add a footer to the app
  tags$footer("© 2023 R Package IntegraTRN. Developed and maintained by Jielin
              Yang."),
  tags$footer("To view the source code or report issues, visit the ",
              tags$a(href = "https://github.com/j-y26/IntegraTRN", 
              "package Github page"), "."),
)


# Define the server for the application
# The server code heavily implements dynamic UI, in which most UI elements are
# generated based on user inputs within the server code
# Note that the same section partitioning is used in the server code as in the
# UI code, which begins with input loading followed by output generation
server <- function(input, output) {

  # === Section 1: select omics data ===========================================
  # RNAseq (required)
  # Small RNAseq (optional)
  # Proteomics (optional)
  # ATACseq (optional)
  # External data
  # miRNA - target gene interactions (optional)
  # TF - target gene interactions (optional)

  # Define a vector of user-selected omics data types
  omicsDataTypes <- reactive({
    # RNAseq is required
    omicsDataTypes <- c(RNA)
    # Add other data types if selected
    if (input$omicSmallRNA) {
      omicsDataTypes <- c(omicsDataTypes, SMALLRNA)
    }
    if (input$omicProteomics) {
      omicsDataTypes <- c(omicsDataTypes, PROTEIN)
    }
    if (input$omicATAC) {
      omicsDataTypes <- c(omicsDataTypes, ATAC)
    }
    if (input$externalmiRNA) {
      omicsDataTypes <- c(omicsDataTypes, MIRNA)
    }
    if (input$externalTF) {
      omicsDataTypes <- c(omicsDataTypes, TF)
    }
    return(omicsDataTypes)
  })

  # Processing data for downloads
  # RNAseq, small RNAseq, and proteomics data are serialized as R objects

  # RNAseq Count Matrix
  output$downloadRnaseqCountMatrix <- downloadHandler(
    filename = function() {
      paste("rnaseq_count_matrix_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(RNAseq_heart, file, row.names = TRUE)
    }
  )
  # RNAseq Sample Metadata
  output$downloadRnaseqSampleMetadata <- downloadHandler(
    filename = function() {
      paste("rnaseq_sample_metadata_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(RNAseq_heart_samples, file, row.names = TRUE)
    }
  )

  # Small RNAseq Count Matrix
  output$downloadSmallRnaseqCountMatrix <- downloadHandler(
    filename = function() {
      paste("small_rnaseq_count_matrix_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(smallRNAseq_heart, file, row.names = TRUE)
    }
  )
  # Small RNAseq Sample Metadata
  output$downloadSmallRnaseqSampleMetadata <- downloadHandler(
    filename = function() {
      paste("small_rnaseq_sample_metadata_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(smallRNAseq_heart_samples, file, row.names = TRUE)
    }
  )

  # Proteomics Count Matrix
  output$downloadProteomicsCountMatrix <- downloadHandler(
    filename = function() {
      paste("proteomics_count_matrix_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(protein_heart, file, row.names = TRUE)
    }
  )
  # Proteomics Sample Metadata
  output$downloadProteomicsSampleMetadata <- downloadHandler(
    filename = function() {
      paste("proteomics_sample_metadata_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(protein_heart_samples, file, row.names = TRUE)
    }
  )
  
  # ATACseq Peaks for Condition 1
  output$downloadATACPeak1 <- downloadHandler(
    filename = function() {
      paste("atac_peaks_condition_1", "bed", sep = ".")
    },
    content = function(file) {
      file.copy(system.file("extdata", "peak1.bed", package = "IntegraTRN"),
                file)
    }
  )
  # ATACseq Peaks for Condition 2
  output$downloadATACPeak2 <- downloadHandler(
    filename = function() {
      paste("atac_peaks_condition_2", "bed", sep = ".")
    },
    content = function(file) {
      file.copy(system.file("extdata", "peak2.bed", package = "IntegraTRN"),
                file)
    }
  )

  # miRNA - Target Interactions
  output$downloadmiRNATarget <- downloadHandler(
    filename = function() {
      paste("mirna_target_interactions", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(as.data.frame(miR2Genes), file, row.names = FALSE)
    }
  )
  
  # TF - Target Interactions
  output$downloadTFTarget <- downloadHandler(
    filename = function() {
      paste("tf_target_interactions", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(as.data.frame(tf2Genes), file, row.names = FALSE)
    }
  )

  # Small RNA Annotation
  output$downloadsmallRNAAnnotation <- downloadHandler(
    filename = function() {
      paste("small_rna_annotation", "csv", sep = ".")
    },
    content = function(file) {
      SNCANNOLIST_HSAPIENS %>%
        utils::stack() %>%
        dplyr::rename(transcript = values, type = ind) %>%
        write.csv(file, row.names = FALSE)
    }
  )

  ##############################################################################
  # Protein name conversion to be implemented later
  ##############################################################################





  # === Section 2: raw data input ==============================================






  
  # RNAseq data input

  # RNAseq count matrix
  rnaseqCountMatrix <- reactive({
    req(input$rnaseqCountMatrix)  # required input

    # Validate file extension and read the file
    ext <- tools::file_ext(input$rnaseqCountMatrix$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    rnaseq <- read.csv(input$rnaseqCountMatrix$datapath, header = TRUE,
                       row.names = 1, stringsAsFactors = FALSE)
    rnaseq <- as.matrix(rnaseq)
    return(rnaseq)
  })

  # RNAseq sample metadata
  rnaseqSampleMetadata <- reactive({
    req(input$rnaseqSampleMetadata)  # required input
    
    # Validate file extension and read the file
    ext <- tools::file_ext(input$rnaseqSampleMetadata$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    rnaseqSampleMetadata <- read.csv(input$rnaseqSampleMetadata$datapath,
                                     header = TRUE, row.names = 1,
                                     stringsAsFactors = FALSE)
    return(rnaseqSampleMetadata)
  })
  # Retrieve the attribute (column) names of the sample metadata
  rnaseqSampleAttributes <- reactive({
    req(input$rnaseqSampleMetadata)
    colnames(rnaseqSampleMetadata())
  })
  # Render selector UI for sample attributes for RNAseq
  output$rnaseqSampleGrouping <- renderUI({
    selectInput(inputId = "rnaseqSampleGrouping",
                label = "Select the grouping variable for RNAseq samples:",
                choices = rnaseqSampleAttributes(),
                multiple = FALSE)
  })

  # Small RNAseq data input
  
  # Small RNAseq count matrix
  smallRnaseqCountMatrix <- reactive({
    if (!is.null(input$smallRnaseqCountMatrix)) { # optional input
      # Validate file extension and read the file
      ext <- tools::file_ext(input$smallRnaseqCountMatrix$datapath)
      validate(need(ext == "csv", "Please upload a .csv file"))
      smallRnaseq <- read.csv(input$smallRnaseqCountMatrix$datapath,
                              header = TRUE, row.names = 1,
                              stringsAsFactors = FALSE)
      smallRnaseq <- as.matrix(smallRnaseq)
      return(smallRnaseq)
    } else {
      return(NULL)
    }
  })

  # Small RNAseq sample metadata
  smallRnaseqSampleMetadata <- reactive({
    if (!is.null(input$smallRnaseqSampleMetadata)) { # optional input
      # Validate file extension and read the file
      ext <- tools::file_ext(input$smallRnaseqSampleMetadata$datapath)
      validate(need(ext == "csv", "Please upload a .csv file"))
      smallRnaseqSampleMetadata <- read.csv(input$smallRnaseqSampleMetadata$datapath,
                                            header = TRUE, row.names = 1,
                                            stringsAsFactors = FALSE)
      return(smallRnaseqSampleMetadata)
    } else {
      return(NULL)
    }
  })

  # Retrieve the attribute (column) names of the sample metadata
  smallRnaseqSampleAttributes <- reactive({
    req(input$smallRnaseqSampleMetadata)
    colnames(smallRnaseqSampleMetadata())
  })
  # Render selector UI for sample attributes for small RNAseq
  output$smallRnaSampleGrouping <- renderUI({
    selectInput(inputId = "smallRnaSampleGrouping",
                label = "Select the grouping variable for small RNA samples:",
                choices = smallRnaseqSampleAttributes(),
                multiple = FALSE)
  })

  # Proteomics data input

  # Proteomics count matrix
  proteomicsCountMatrix <- reactive({
    if (!is.null(input$proteomicsCountMatrix)) { # optional input
      # Validate file extension and read the file
      ext <- tools::file_ext(input$proteomicsCountMatrix$datapath)
      validate(need(ext == "csv", "Please upload a .csv file"))
      proteomics <- read.csv(input$proteomicsCountMatrix$datapath,
                             header = TRUE, row.names = 1,
                             stringsAsFactors = FALSE)
      proteomics <- as.matrix(proteomics)
      return(proteomics)
    } else {
      return(NULL)
    }
  })

  # Proteomics sample metadata
  proteomicsSampleMetadata <- reactive({
    if (!is.null(input$proteomicsSampleMetadata)) { # optional input
      # Validate file extension and read the file
      ext <- tools::file_ext(input$proteomicsSampleMetadata$datapath)
      validate(need(ext == "csv", "Please upload a .csv file"))
      proteomicsSampleMetadata <- read.csv(input$proteomicsSampleMetadata$datapath,
                                           header = TRUE, row.names = 1,
                                           stringsAsFactors = FALSE)
      return(proteomicsSampleMetadata)
    } else {
      return(NULL)
    }
  })
  # Retrieve the attribute (column) names of the sample metadata
  proteomicsSampleAttributes <- reactive({
    req(input$proteomicsSampleMetadata)
    colnames(proteomicsSampleMetadata())
  })
  # Render selector UI for sample attributes for proteomics
  output$proteomicsSampleGrouping <- renderUI({
    selectInput(inputId = "proteomicsSampleGrouping",
                label = "Select the grouping variable for proteomics samples:",
                choices = proteomicsSampleAttributes(),
                multiple = FALSE)
  })



  





  # Integrative input validation
  # Check that all provided inputs satisfy the requirements for integrative
  # analysis, in addition to the specific requirements for each data type







  # ============================================================================
  # === FOR OUTPUT TESTING ONLY, REMOVE LATER ==================================
  # ============================================================================

  # Check that the omics data types are correctly selected
  output$omicsDataTypes <- renderText({
    paste(omicsDataTypes(), collapse = ", ")
  })






  
}



# Build the shiny app object
shinyApp(ui = ui, server = server)

# [END]