# Purpose: Define the shiny app for the package
# Author: Jielin Yang
# Date: 2023-11-29
# Version: 1.0
# Bugs and Issues: None

# Load libraries
library(shiny)
library(shinyalert)
library(shinyBS)

# Define some global variables
# Note that RNA, SMALLRNA, ATAC, and PROTEIN is already defined within the
# package
# The following variables are used to define external data types
# Note these definitions are only available within the app
MIRNA <- "miRNA_target"
TF <- "TF_target"
HSAPIENS <- "Human (hg38)"
MMUSCULUS <- "Mouse (mm10)"


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
    the inference of core transcriptional regulatory networks")),

  # Use sidebar layout with a sidebar panel and a main panel
  sidebarLayout(

    # Sidebar panel for inputs
    sidebarPanel(

      # Creating multiple tabs for different input sections
      tabsetPanel(

        # === Panel for introduction ===========================================
        tabPanel(
          "Introduction",
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
          tags$img(
            src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/Schematics.jpg", # link cannot be broken to shorten the line
            width = "550px", align = "center"
          ),
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
        tabPanel(
          "Section 1",
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
          column(
            width = 6,
            tags$img(
              src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_rnaseq.jpg",
              width = "200px", align = "center"
            ),
            div(style = "margin-top: 10px"),
            tags$b(id = "omicRNA", "RNAseq (required)"), align = "center"
          ),
          bsTooltip(
            id = "omicRNA",
            title = paste0(
              "Bulk RNA sequencing (RNAseq) data is ",
              "required for the analysis, which provides a full ",
              "transcriptomic profile of gene expression for the ",
              "tissue/biological system of interest."
            )
          ),
          # The use of bsTooltip is to add hover-over tooltips is derived from
          # https://stackoverflow.com/questions/16449252/tooltip-on-shiny-r
          # by Praneeth Jagarapu

          # Small RNAseq
          column(
            width = 6,
            tags$img(
              src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_smallrnaseq.jpg",
              width = "200px", align = "center"
            ),
            checkboxInput(
              inputId = "omicSmallRNA",
              label = tags$b("Small RNAseq"),
              value = FALSE
            ), align = "center"
          ),
          bsTooltip(
            id = "omicSmallRNA",
            title = paste0(
              "Small RNA sequencing (small RNAseq) ",
              "data is optional for the analysis, which provides ",
              "insights into the expression of small RNAs, such ",
              "as miRNAs, piRNAs, and snoRNAs. They are important ",
              "regulators of gene expression."
            )
          ),

          # Proteomics
          column(
            width = 6,
            tags$img(
              src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_proteomics.jpg",
              width = "200px", align = "center"
            ),
            checkboxInput(
              inputId = "omicProteomics",
              label = tags$b("Proteomics"),
              value = FALSE
            ), align = "center"
          ),
          bsTooltip(
            id = "omicProteomics",
            title = paste0(
              "Proteomics data is optional for the ",
              "analysis. Cross referencing transcriptomic and ",
              "proteomic data can further validate gene ",
              "expression changes at the protein level."
            )
          ),

          # ATACseq
          column(
            width = 6,
            tags$img(
              src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_atac.jpg",
              width = "200px", align = "center"
            ),
            checkboxInput(
              inputId = "omicATAC",
              label = tags$b("ATACseq"),
              value = FALSE
            ), align = "center"
          ),
          bsTooltip(
            id = "omicATAC",
            title = paste0(
              "ATACseq data is optional. It ",
              "provides valuable information on chromatin ",
              "accessibility, which reflects putative ",
              "binding sites of transcription factors. "
            )
          ),

          # External data
          # miRNA - target gene interactions
          column(
            width = 6,
            tags$img(
              src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_mirna.jpg",
              width = "200px", align = "center"
            ),
            checkboxInput(
              inputId = "externalmiRNA",
              label = tags$b("microRNA - Target Interactions"),
              value = FALSE
            ), align = "center"
          ),
          bsTooltip(
            id = "externalmiRNA",
            title = paste0(
              "microRNAs (miRNAs) are important ",
              "regulators of gene expression by binding to ",
              "target genes and leading to their degradation ",
              "or translational repression. Curated miRNA - ",
              "target gene interactions are available from ",
              "external databases, and can be used to establish ",
              "a non-tissue-specific context for the analysis."
            )
          ),

          # TF - target gene interactions
          column(
            width = 6,
            tags$img(
              src = "https://raw.githubusercontent.com/j-y26/IntegraTRN/master/inst/extdata/icon_tf.jpg",
              width = "200px", align = "center"
            ),
            checkboxInput(
              inputId = "externalTF",
              label = tags$b("TF - Target Interactions"),
              value = FALSE
            ), align = "center"
          ),
          bsTooltip(
            id = "externalTF",
            title = paste0(
              "Transcription factors (TFs) ",
              "can mediate transcriptional activation or ",
              "repression of target genes. Curated TF - ",
              "target gene interactions are available from ",
              "previous ChIP-seq studies, but may lack ",
              "tissue-specificity. Intersecting with ATACseq ",
              "data can help to identify tissue-specific ",
              "TF - target gene interactions."
            )
          ),
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
            tags$li(downloadLink(
              outputId = "downloadRnaseqCountMatrix",
              label = "RNAseq count matrix (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadRnaseqSampleMetadata",
              label = "RNAseq sample metadata (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadSmallRnaseqCountMatrix",
              label = "Small RNAseq count matrix (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadSmallRnaseqSampleMetadata",
              label = "Small RNAseq sample metadata (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadProteomicsCountMatrix",
              label = "Proteomics count matrix (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadProteomicsSampleMetadata",
              label = "Proteomics sample metadata (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadATACPeak1",
              label = "ATACseq peak for condition 1 (.bed)"
            )),
            tags$li(downloadLink(
              outputId = "downloadATACPeak2",
              label = "ATACseq peak for condition 2 (.bed)"
            )),
            tags$li(downloadLink(
              outputId = "downloadmiRNATarget",
              label = "miRNA - target interactions (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadTFTarget",
              label = "TF - target interactions (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadsmallRNAAnnotation",
              label = "Small RNA annotation (.csv)"
            )),
            tags$li(downloadLink(
              outputId = "downloadproteinNameConversion",
              label = "Gene - protein name conversion (.csv)"
            ))
          ),
        ),


        # === Panel for section 2: raw data input ==============================
        tabPanel(
          "Section 2",
          tags$h4("Section 2: Raw Data Input"),
          br(),

          # Description of the section
          tags$p(tags$b("Description:"), "In this section, please upload
          the raw data files for each data type according to the selected
          data types in the previous section."),
          br(),

          # RNAseq data input (always shown)
          tags$h5("Part 1: RNAseq Data Input (Required)"),

          # RNAseq count matrix
          tags$b("Upload RNAseq count matrix (.csv):"),
          fileInput(
            inputId = "rnaseqCountMatrix",
            label = "A numeric matrix with genes as rows
                  and samples as columns. The gene names should be specified in
                  the first column. The first row should be sample names.",
            accept = c(".csv")
          ),
          div(style = "margin-top: -35px"), # move the download link up
          # The above solution to adjust the position of the download link is
          # from https://stackoverflow.com/questions/62476054/reduce-space-
          # between-fileinput-and-text-in-shiny by Stéphane Laurent
          downloadLink(
            outputId = "onsiteDownloadRnaseqCountMatrix",
            label = "Download example RNAseq count matrix (.csv)"
          ),
          br(),

          # RNAseq sample metadata
          tags$b("Upload RNAseq sample metadata (.csv):"),
          fileInput(
            inputId = "rnaseqSampleMetadata",
            label = "A data frame with samples as rows and attributes
                  as columns. The first column should be sample names. The
                  first row should be metadata names.",
            accept = c(".csv")
          ),
          div(style = "margin-top: -35px"), # move the download link up
          downloadLink(
            outputId = "onsiteDownloadRnaseqSampleMetadata",
            label = "Download example RNAseq sample metadata (.csv)"
          ),
          br(),
          # Select the grouping variable for RNAseq samples
          uiOutput("rnaseqSampleGrouping"),
          br(),
          br(),

          # Small RNAseq data input (optional)
          uiOutput("smallRnaRawData"),

          # Proteomics data input (optional)
          uiOutput("proteomicsRawData"),

          # ATACseq data input (optional)
          uiOutput("atacRawData"),
        ),

        # === Panel for section 3: annotations and external data ===============
        tabPanel(
          "Section 3",
          tags$h4("Section 3: Annotations and Externally Curated Regulatory
                Interactions"),
          br(),

          # Description of the section
          tags$p(tags$b("Description:"), "In this section, please upload
          annotation files for small RNAs and proteins, as well as
          externally curated regulatory interactions, such as TF-target
          gene and miRNA-target gene interactions."),
          tags$p(
            "An example tool to search for externally curated interactions is
                miRNet, available ",
            tags$a(href = "https://www.mirnet.ca/", "miRNet"), "."
          ),
          br(),
          br(),
          tags$h5("Part 1: Annotations"),
          selectInput(
            inputId = "genAssembly",
            label = "Please select the genome assembly:",
            choices = c(HSAPIENS, MMUSCULUS),
            selected = HSAPIENS
          ),

          # Small RNA annotation
          uiOutput("smallRNAAnnotation"),

          # Protein name conversion
          uiOutput("proteinNameConversion"),

          # Check annotation coverage button
          uiOutput("checkAnnotationCoverage"),
          bsAlert(anchorId = "smallRNACoverageAlert"),
          bsAlert(anchorId = "smallRNAFullCoverageAlert"),
          bsAlert(anchorId = "proteinCoverageAlert"),
          bsAlert(anchorId = "proteinFullCoverageAlert"),
          br(),

          # External data
          uiOutput("externalRawData"),

          # miRNA - target interactions
          uiOutput("miRNATargetRawData"),

          # TF - target interactions
          uiOutput("tfTargetRawData"),





          # Space holder, implement later
        ),

        #   # === Panel for section 4: integrative analysis ====================
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
        tabPanel(
          "Output Testing",
          tags$h4("Output Testing"),
          br(),

          # Test correct selection of omics data types
          tags$p("The following text should be the selected omics data types:"),
          textOutput(outputId = "omicsDataTypes"),
          br(),

          # Test correct input of RNAseq data
          tags$p("RNAseq raw data input:"),
          tableOutput(outputId = "rnaCM"),
          tableOutput(outputId = "rnaSM"),
          textOutput(outputId = "rnaSG"),

          # Test correct input of small RNAseq data
          tags$p("Small RNAseq raw data input:"),
          tableOutput(outputId = "smallRnaCM"),
          tableOutput(outputId = "smallRnaSM"),
          textOutput(outputId = "smallRnaSG"),

          # Test correct input of proteomics data
          tags$p("Proteomics raw data input:"),
          tableOutput(outputId = "proteomicsCM"),
          tableOutput(outputId = "proteomicsSM"),
          textOutput(outputId = "proteomicsSG"),

          # Test correct input of ATACseq data
          tags$p("ATACseq raw data input:"),
          tableOutput(outputId = "atacPeak1Test"),
          tableOutput(outputId = "atacPeak2Test"),

          # Test correct input of miRNA - target interactions
          tags$p("miRNA - target interactions:"),
          tableOutput(outputId = "miRNATargetTest"),

          # Test correct input of TF - target interactions
          tags$p("TF - target interactions:"),
          tableOutput(outputId = "TFTargetTest"),

          # Test correct input of small RNA annotation
          tags$p("Small RNA annotation:"),
          tableOutput(outputId = "smallRNAAnnotationTest"),

          # Test correct input of protein name conversion
          tags$p("Protein name conversion:"),
          tableOutput(outputId = "proteinNameConversionTest"),

          # Test correct genome assembly selection
          tags$p("Genome assembly:"),
          textOutput(outputId = "genAssemblyTest"),
        ),
      ),
    ),
  ),

  # Add a footer to the app
  tags$footer("© 2023 R Package IntegraTRN. Developed and maintained by Jielin
              Yang."),
  tags$footer(
    "To view the source code or report issues, visit the ",
    tags$a(
      href = "https://github.com/j-y26/IntegraTRN",
      "package Github page"
    ), "."
  ),
)


# Define the server for the application
# The server code heavily implements dynamic UI, in which most UI elements are
# generated based on user inputs within the server code
# Note that the same section partitioning is used in the server code as in the
# UI code, which begins with input loading followed by output generation
server <- function(input, output, session) {
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
  downloadRnaseqCountMatrix <- downloadHandler(
    filename = function() {
      paste("rnaseq_count_matrix_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(RNAseq_heart, file, row.names = TRUE)
    }
  )
  output$downloadRnaseqCountMatrix <- downloadRnaseqCountMatrix

  # RNAseq Sample Metadata
  downloadRnaseqSampleMetadata <- downloadHandler(
    filename = function() {
      paste("rnaseq_sample_metadata_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(RNAseq_heart_samples, file, row.names = TRUE)
    }
  )
  output$downloadRnaseqSampleMetadata <- downloadRnaseqSampleMetadata

  # Small RNAseq Count Matrix
  downloadSmallRnaseqCountMatrix <- downloadHandler(
    filename = function() {
      paste("small_rnaseq_count_matrix_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(smallRNAseq_heart, file, row.names = TRUE)
    }
  )
  output$downloadSmallRnaseqCountMatrix <- downloadSmallRnaseqCountMatrix

  # Small RNAseq Sample Metadata
  downloadSmallRnaseqSampleMetadata <- downloadHandler(
    filename = function() {
      paste("small_rnaseq_sample_metadata_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(smallRNAseq_heart_samples, file, row.names = TRUE)
    }
  )
  output$downloadSmallRnaseqSampleMetadata <- downloadSmallRnaseqSampleMetadata

  # Proteomics Count Matrix
  downloadProteomicsCountMatrix <- downloadHandler(
    filename = function() {
      paste("proteomics_count_matrix_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(protein_heart, file, row.names = TRUE)
    }
  )
  output$downloadProteomicsCountMatrix <- downloadProteomicsCountMatrix

  # Proteomics Sample Metadata
  downloadProteomicsSampleMetadata <- downloadHandler(
    filename = function() {
      paste("proteomics_sample_metadata_fetal_heart", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(protein_heart_samples, file, row.names = TRUE)
    }
  )
  output$downloadProteomicsSampleMetadata <- downloadProteomicsSampleMetadata

  # ATACseq Peaks for Condition 1
  downloadATACPeak1 <- downloadHandler(
    filename = function() {
      paste("atac_peaks_condition_1", "bed", sep = ".")
    },
    content = function(file) {
      file.copy(
        system.file("extdata", "peak1.bed", package = "IntegraTRN"),
        file
      )
    }
  )
  output$downloadATACPeak1 <- downloadATACPeak1

  # ATACseq Peaks for Condition 2
  downloadATACPeak2 <- downloadHandler(
    filename = function() {
      paste("atac_peaks_condition_2", "bed", sep = ".")
    },
    content = function(file) {
      file.copy(
        system.file("extdata", "peak2.bed", package = "IntegraTRN"),
        file
      )
    }
  )
  output$downloadATACPeak2 <- downloadATACPeak2

  # miRNA - Target Interactions
  downloadmiRNATarget <- downloadHandler(
    filename = function() {
      paste("mirna_target_interactions", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(as.data.frame(miR2Genes), file, row.names = FALSE)
    }
  )
  output$downloadmiRNATarget <- downloadmiRNATarget

  # TF - Target Interactions
  downloadTFTarget <- downloadHandler(
    filename = function() {
      paste("tf_target_interactions", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(as.data.frame(tf2Genes), file, row.names = FALSE)
    }
  )
  output$downloadTFTarget <- downloadTFTarget

  # Small RNA Annotation
  downloadsmallRNAAnnotation <- downloadHandler(
    filename = function() {
      paste("small_rna_annotation_hsapiens", "csv", sep = ".")
    },
    content = function(file) {
      SNCANNOLIST_HSAPIENS %>%
        utils::stack() %>%
        dplyr::rename(transcript = values, type = ind) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  output$downloadsmallRNAAnnotation <- downloadsmallRNAAnnotation

  # Protein to Gene Name Conversion
  downloadproteinNameConversion <- downloadHandler(
    filename = function() {
      paste("protein_name_conversion_partial", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(proteinGeneIDConvert, file, row.names = FALSE)
    }
  )
  output$downloadproteinNameConversion <- downloadproteinNameConversion





  # === Section 2: raw data input ==============================================

  # RNAseq data (required)
  # RNAseq count matrix
  rnaseqCountMatrix <- reactive({
    req(input$rnaseqCountMatrix)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$rnaseqCountMatrix$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    rnaseq <- read.csv(input$rnaseqCountMatrix$datapath,
      header = TRUE,
      row.names = 1, stringsAsFactors = FALSE
    )
    rnaseq <- as.matrix(rnaseq)
    return(rnaseq)
  })

  # Processing onsite download option
  output$onsiteDownloadRnaseqCountMatrix <- downloadRnaseqCountMatrix

  # RNAseq sample metadata
  rnaseqSampleMetadata <- reactive({
    req(input$rnaseqSampleMetadata)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$rnaseqSampleMetadata$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    rnaseqSampleMetadata <- read.csv(input$rnaseqSampleMetadata$datapath,
      header = TRUE, row.names = 1,
      stringsAsFactors = FALSE
    )
    return(rnaseqSampleMetadata)
  })
  # Processing onsite download option
  output$onsiteDownloadRnaseqSampleMetadata <- downloadRnaseqSampleMetadata

  # Retrieve the attribute (column) names of the sample metadata
  rnaseqSampleAttributes <- reactive({
    req(input$rnaseqSampleMetadata)
    colnames(rnaseqSampleMetadata())
  })
  # Render selector UI for sample attributes for RNAseq
  output$rnaseqSampleGrouping <- renderUI({
    selectInput(
      inputId = "rnaseqSampleGrouping",
      label = "Select the grouping variable for RNAseq samples:",
      choices = rnaseqSampleAttributes(),
      multiple = FALSE
    )
  })

  # Small RNAseq data (dynamic UI)
  useSmallRNAseq <- reactive({
    return(SMALLRNA %in% omicsDataTypes())
  })

  output$smallRnaRawData <- renderUI({
    if (useSmallRNAseq()) {
      # Calculate a part number
      partNumber <- which(omicsDataTypes() == SMALLRNA)

      tagList(
        # Small RNAseq count matrix
        tags$h5("Part ", partNumber, ": Small RNAseq Data Input"),
        tags$b("Upload small RNAseq count matrix (.csv):"),
        fileInput(
          inputId = "smallRnaseqCountMatrix",
          label = "A numeric matrix with transcripts as rows
                and samples as columns. The small RNA names should be specified
                in the first column. The first row should be sample names.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadSmallRnaseqCountMatrix",
          label = "Download example small RNAseq count matrix (.csv)"
        ),
        br(),
        br(),

        # Small RNAseq sample metadata
        tags$b("Upload small RNAseq sample metadata (.csv):"),
        fileInput(
          inputId = "smallRnaseqSampleMetadata",
          label = "A data frame with samples as rows and attributes
                  as columns. The first column should be sample names. The
                  first row should be metadata names.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadSmallRnaseqSampleMetadata",
          label = "Download example small RNAseq sample metadata (.csv)"
        ),
        br(),
        # Select the grouping variable for small RNAseq samples
        uiOutput("smallRnaSampleGrouping"),
        br(),
        br(),
      )
    }
  })

  # Small RNAseq count matrix (server)
  smallRnaseqCountMatrix <- reactive({
    req(input$smallRnaseqCountMatrix)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$smallRnaseqCountMatrix$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    smallRnaseq <- read.csv(input$smallRnaseqCountMatrix$datapath,
      header = TRUE,
      row.names = 1, stringsAsFactors = FALSE
    )
    smallRnaseq <- as.matrix(smallRnaseq)
    return(smallRnaseq)
  })

  # Processing onsite download option
  output$onsiteDownloadSmallRnaseqCountMatrix <- downloadSmallRnaseqCountMatrix

  # Small RNAseq sample metadata (server)
  smallRnaseqSampleMetadata <- reactive({
    req(input$smallRnaseqSampleMetadata)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$smallRnaseqSampleMetadata$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    smallRnaseqSampleMetadata <- read.csv(input$smallRnaseqSampleMetadata$datapath,
      header = TRUE, row.names = 1,
      stringsAsFactors = FALSE
    )
    return(smallRnaseqSampleMetadata)
  })

  # Processing onsite download option
  output$onsiteDownloadSmallRnaseqSampleMetadata <-
    downloadSmallRnaseqSampleMetadata

  # Retrieve the attribute (column) names of the sample metadata
  smallRnaseqSampleAttributes <- reactive({
    req(input$smallRnaseqSampleMetadata)
    colnames(smallRnaseqSampleMetadata())
  })
  # Render selector UI for sample attributes for small RNAseq
  output$smallRnaSampleGrouping <- renderUI({
    selectInput(
      inputId = "smallRnaSampleGrouping",
      label = "Select the grouping variable for small RNAseq samples:",
      choices = smallRnaseqSampleAttributes(),
      multiple = FALSE
    )
  })

  # Proteomics data (dynamic UI)
  useProteomics <- reactive({
    return(PROTEIN %in% omicsDataTypes())
  })

  output$proteomicsRawData <- renderUI({
    if (useProteomics()) {
      # Calculate a part number
      partNumber <- which(omicsDataTypes() == PROTEIN)

      tagList(
        # Proteomics count matrix
        tags$h5("Part ", partNumber, ": Proteomics Data Input"),
        tags$b("Upload proteomics count matrix (.csv):"),
        fileInput(
          inputId = "proteomicsCountMatrix",
          label = "A numeric matrix with proteins as rows
                  and samples as columns. The protein names should be specified
                  in the first column. The first row should be sample names.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadProteomicsCountMatrix",
          label = "Download example proteomics count matrix (.csv)"
        ),
        br(),
        br(),

        # Proteomics sample metadata
        tags$b("Upload proteomics sample metadata (.csv):"),
        fileInput(
          inputId = "proteomicsSampleMetadata",
          label = "A data frame with samples as rows and attributes
                  as columns. The first column should be sample names. The
                  first row should be metadata names.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadProteomicsSampleMetadata",
          label = "Download example proteomics sample metadata (.csv)"
        ),
        br(),
        # Select the grouping variable for proteomics samples
        uiOutput("proteomicsSampleGrouping"),
        br(),
        br(),
      )
    }
  })

  # Proteomics count matrix (server)
  proteomicsCountMatrix <- reactive({
    req(input$proteomicsCountMatrix)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$proteomicsCountMatrix$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    proteomics <- read.csv(input$proteomicsCountMatrix$datapath,
      header = TRUE,
      row.names = 1, stringsAsFactors = FALSE
    )
    proteomics <- as.matrix(proteomics)
    return(proteomics)
  })

  # Processing onsite download option
  output$onsiteDownloadProteomicsCountMatrix <- downloadProteomicsCountMatrix

  # Proteomics sample metadata (server)
  proteomicsSampleMetadata <- reactive({
    req(input$proteomicsSampleMetadata)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$proteomicsSampleMetadata$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    proteomicsSampleMetadata <- read.csv(input$proteomicsSampleMetadata$datapath,
      header = TRUE, row.names = 1,
      stringsAsFactors = FALSE
    )
    return(proteomicsSampleMetadata)
  })

  # Processing onsite download option
  output$onsiteDownloadProteomicsSampleMetadata <-
    downloadProteomicsSampleMetadata

  # Retrieve the attribute (column) names of the sample metadata
  proteomicsSampleAttributes <- reactive({
    req(input$proteomicsSampleMetadata)
    colnames(proteomicsSampleMetadata())
  })
  # Render selector UI for sample attributes for proteomics
  output$proteomicsSampleGrouping <- renderUI({
    selectInput(
      inputId = "proteomicsSampleGrouping",
      label = "Select the grouping variable for proteomics samples:",
      choices = proteomicsSampleAttributes(),
      multiple = FALSE
    )
  })

  # ATACseq data (dynamic UI)
  useATACseq <- reactive({
    return(ATAC %in% omicsDataTypes())
  })

  # Define a set of file formats to be accepted
  peakFormats <- c(".bed", ".narrowPeak", ".broadPeak", ".gappedPeak")

  output$atacRawData <- renderUI({
    if (useATACseq()) {
      # Calculate a part number
      partNumber <- which(omicsDataTypes() == ATAC)

      tagList(
        # ATACseq peaks for condition 1
        tags$h5("Part ", partNumber, ": ATACseq Data Input)"),
        tags$b("Upload ATACseq peaks for condition 1 (.bed):"),
        fileInput(
          inputId = "atacPeak1",
          label = "A .bed file with peaks as rows and samples as columns.
                The first column should be peak names. The first row should be
                sample names.",
          accept = peakFormats
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadATACPeak1",
          label = "Download example ATACseq peaks for condition 1 (.bed)"
        ),
        br(),
        br(),

        # ATACseq peaks for condition 2
        tags$b("Upload ATACseq peaks for condition 2 (.bed):"),
        fileInput(
          inputId = "atacPeak2",
          label = "A .bed file with peaks as rows and samples as columns.
                The first column should be peak names. The first row should be
                sample names.",
          accept = peakFormats
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadATACPeak2",
          label = "Download example ATACseq peaks for condition 2 (.bed)"
        ),
        br(),
        br(),
      )
    }
  })

  # ATACseq peaks for condition 1 (server)
  atacPeak1 <- reactive({
    req(input$atacPeak1)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$atacPeak1$datapath)
    validate(need(ext == "bed", "Please upload a .bed file"))
    atacPeak1 <- read.table(input$atacPeak1$datapath, header = FALSE)
    atacPeak1 <- atacPeak1[, 1:3] %>% setNames(CHROMINFO)
    return(atacPeak1) # keep as data frame
  })

  # Processing onsite download option
  output$onsiteDownloadATACPeak1 <- downloadATACPeak1

  # ATACseq peaks for condition 2 (server)
  atacPeak2 <- reactive({
    req(input$atacPeak2)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$atacPeak2$datapath)
    validate(need(ext == "bed", "Please upload a .bed file"))
    atacPeak2 <- read.table(input$atacPeak2$datapath, header = FALSE)
    atacPeak2 <- atacPeak2[, 1:3] %>% setNames(CHROMINFO)
    return(atacPeak2) # keep as data frame
  })

  # Processing onsite download option
  output$onsiteDownloadATACPeak2 <- downloadATACPeak2




  # === Section 3: annotations and external data ===============================

  # Part 1: annotations

  # Selecting the genome assembly
  genAssembly <- reactive({
    req(input$genAssembly)
    return(input$genAssembly)
  })

  # Small RNA annotation (dynamic UI)
  output$smallRNAAnnotation <- renderUI({
    if (useSmallRNAseq()) {
      tagList(
        # User should upload a small RNA annotation file
        tags$b("Upload small RNA annotation (.csv):"),
        fileInput(
          inputId = "smallRNAAnnotation",
          label = "A .csv file with the first column as small RNA names
                and the second column as the type of small RNAs. Valid types
                include 'miRNA', 'piRNA', 'tRNA', 'circRNA', 'snRNA',
                and 'snoRNA'. The first row should be column names: 'transcript'
                and 'type'.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -20px"),
        tags$b(id = "sncAnnoNote", "Note on small RNA annotation
                                    (click here to see details):"),
        bsPopover(
          id = "sncAnnoNote",
          title = "Human Small RNA Annotations",
          content = paste0(
            "The IntegraTRN package provides a pre-compiled ",
            "small RNA annotation file for human following the ",
            "standard hsa naming convention. This annotation ",
            "utilizes data from miRbase, piRNAbank, piRBase, ",
            "GtRNAdb, circBase, GENCODE, and piRNACluster. The ",
            "annotation file is compiled using the ",
            tags$a(
              href = "https://github.com/cougarlj/COMPSRA",
              "COMPSRA tool"
            ), "."
          ),
          placement = "right",
          trigger = "click"
        ),
        tags$p("A pre-compiled small RNA annotation file for human is available
                for download:"),
        div(style = "margin-top: -10px"),
        downloadLink(
          outputId = "onsiteDownloadsmallRNAAnnotation",
          label = "Download small RNA annotation (.csv)"
        ),
        br(),
        br(),
      )
    }
  })

  # Small RNA annotation (server)
  # First provide a download link
  output$onsiteDownloadsmallRNAAnnotation <- downloadsmallRNAAnnotation

  # Retrieve the small RNA annotation
  smallRNAAnnotation <- reactive({
    req(input$smallRNAAnnotation)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$smallRNAAnnotation$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    smallRNAAnnotation <- read.csv(input$smallRNAAnnotation$datapath,
      header = TRUE, stringsAsFactors = FALSE
    )
    smallRNAAnnotation <- smallRNAAnnotation[, 1:2] %>%
      setNames(c("transcript", "type"))
    return(smallRNAAnnotation)
  })

  # Protein name conversion (dynamic UI)
  output$proteinNameConversion <- renderUI({
    if (useProteomics()) {
      tagList(
        # User should upload a protein name conversion file
        tags$b("Upload gene - protein name conversion (.csv):"),
        fileInput(
          inputId = "proteinNameConversion",
          label = "A .csv file with the first column as protein names
                and the second column as gene names. The first row should
                be column names: 'gene' and 'protein'. Use consistent naming
                with the RNAseq and proteomics data.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -20px"),
        tags$b(
          id = "proteinNameConversionNote",
          "Note on ID conversion (click here to see details):"
        ),
        bsPopover(
          id = "proteinNameConversionNote",
          title = "Converting protein IDs to gene IDs",
          content = paste0(
            "ID mapping is a key issue in bioinformatics. It ",
            "is recommended that the users use official tools that ",
            "provide a comprehensive mapping between gene and ",
            "protein IDs. It is required that the IDs of genes and ",
            "proteins are consistent with the RNAseq and proteomics ",
            "data. The provided example file illustrates an ",
            "non-comprehensive case, where the example does not ",
            "fully cover the proteins in the proteomics data. The ",
            "programs can handle these cases by ignoring the ",
            "unmapped proteins, but the constructed TRN may lose ",
            "some information."
          ),
          placement = "right",
          trigger = "click"
        ),
        tags$p("An example gene - protein name conversion file is available
                for download:"),
        div(style = "margin-top: -10px"),
        downloadLink(
          outputId = "onsiteDownloadproteinNameConversion",
          label = "Download gene - protein name conversion (.csv)"
        ),
        br(),
        br(),
      )
    }
  })

  # Protein name conversion (server)
  # First provide a download link
  output$onsiteDownloadproteinNameConversion <- downloadproteinNameConversion

  # Retrieve the protein name conversion
  proteinNameConversion <- reactive({
    req(input$proteinNameConversion)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$proteinNameConversion$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    proteinNameConversion <- read.csv(input$proteinNameConversion$datapath,
      header = TRUE, stringsAsFactors = FALSE
    )
    proteinNameConversion <- proteinNameConversion[, 1:2] %>%
      setNames(c("protein", "gene"))
    return(proteinNameConversion)
  })

  # Button to check for annotation coverage (dynamic UI)
  output$checkAnnotationCoverage <- renderUI({
    if (useSmallRNAseq() || useProteomics()) {
      tagList(
        actionButton(
          inputId = "checkAnnotationCoverage",
          label = "Check Annotation Coverage"
        ),
        br(),
        br(),
      )
    }
  })

  # Check annotation coverage
  # Check if the small RNA annotation covers all small RNAs in the small RNAseq
  # data
  # Check if the protein name conversion covers all proteins in the proteomics
  # data
  # If not, display a warning message
  observeEvent(input$checkAnnotationCoverage, {
    if (useSmallRNAseq()) {
      # Check if the small RNA annotation covers all small RNAs in the small
      # RNAseq data
      smallRNAs <- rownames(smallRnaseqCountMatrix())
      if (!all(smallRNAs %in% smallRNAAnnotation()[, 1])) {
        createAlert(session,
          anchorId = "smallRNACoverageAlert",
          alertId = "smallRNACoverageAlertID",
          title = "Incomplete Small RNA Annotation Coverage",
          content = paste0(
            "The small RNA annotation does not cover all small RNAs in the ",
            "small RNAseq data. Unannotated small RNAs will not be factored ",
            "into the TRN construction."
          )
        )
      } else {
        createAlert(session,
          anchorId = "smallRNAFullCoverageAlert",
          alertId = "smallRNAFullCoverageAlertID",
          title = "Great! All small RNAs are annotated."
        )
      }
    }
    if (useProteomics()) {
      # Same for the protein name conversion
      proteins <- rownames(proteomicsCountMatrix())
      if (!all(proteins %in% proteinNameConversion()$protein)) {
        createAlert(session,
          anchorId = "proteinCoverageAlert",
          alertId = "proteinCoverageAlertID",
          title = "Incomplete Protein Coverage",
          content = paste0(
            "The protein-gene ID conversion infomration does not cover all ", "proteins in the ",
            "proteomics data. Uncovered proteins and their corresoinding ",
            "genes may be lost during TRN construction."
          )
        )
      } else {
        createAlert(session,
          anchorId = "proteinFullCoverageAlert",
          alertId = "proteinFullCoverageAlertID",
          title = "Great! All proteins are covered."
        )
      }
    }
  })



  # Part 2: external interaction data

  # External raw data header
  output$externalRawData <- renderUI({
    if (usemiRNATarget() || useTFTarget()) {
      tagList(
        tags$h5("Part 2: Externally Curated Regulatory Interactions"),
      )
    }
  })


  # External miRNA - target interactions (dynamic UI)
  usemiRNATarget <- reactive({
    return(MIRNA %in% omicsDataTypes())
  })

  output$miRNATargetRawData <- renderUI({
    if (usemiRNATarget()) {
      tagList(
        # miRNA - target interactions
        tags$b("Upload miRNA - target interactions (.csv):"),
        fileInput(
          inputId = "miRNATarget",
          label = "A .csv file with the first column as miRNA names and
                the second column as target gene names. The naming should be
                consistent with the RNAseq and small RNAseq data. The first row
                should be column names: 'regulator' and 'target'.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadmiRNATarget",
          label = "Download example miRNA - target interactions (.csv)"
        ),
        br(),
        br(),
      )
    }
  })

  # miRNA - target interactions (server)
  miRNATarget <- reactive({
    req(input$miRNATarget)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$miRNATarget$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    miRNATarget <- read.csv(input$miRNATarget$datapath,
      header = TRUE,
      stringsAsFactors = FALSE
    )
    miRNATarget <- miRNATarget[, 1:2] %>% setNames(c("regulator", "target"))
    return(miRNATarget)
  })

  # Processing onsite download option
  output$onsiteDownloadmiRNATarget <- downloadmiRNATarget

  # External TF - target interactions (dynamic UI)
  useTFTarget <- reactive({
    return(TF %in% omicsDataTypes())
  })

  output$tfTargetRawData <- renderUI({
    if (useTFTarget()) {
      # Calculate a part number
      partNumber <- which(omicsDataTypes() == TF)
      tagList(
        # TF - target interactions
        tags$b("Upload TF - target interactions (.csv):"),
        fileInput(
          inputId = "tfTarget",
          label = "A .csv file with the first column as transcription factors
            the second column as target gene names. The naming should be
            consistent with the RNAseq data. The first row
            should be column names: 'regulator' and 'target'.",
          accept = c(".csv")
        ),
        div(style = "margin-top: -35px"), # move the download link up
        downloadLink(
          outputId = "onsiteDownloadTFTarget",
          label = "Download example TF - target interactions (.csv)"
        ),
        br(),
        br(),
      )
    }
  })

  # TF - target interactions (server)
  TFTarget <- reactive({
    req(input$tfTarget)

    # Validate file extension and read the file
    ext <- tools::file_ext(input$tfTarget$datapath)
    validate(need(ext == "csv", "Please upload a .csv file"))
    TFTarget <- read.csv(input$tfTarget$datapath,
      header = TRUE,
      stringsAsFactors = FALSE
    )
    TFTarget <- TFTarget[, 1:2] %>% setNames(c("regulator", "target"))
    return(TFTarget)
  })

  # Processing onsite download option
  output$onsiteDownloadTFTarget <- downloadTFTarget


















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

  # Display the first 2 rows of the RNAseq count matrix, sample metadata, and
  # the selected grouping variable
  output$rnaCM <- renderTable({
    rnaseqCountMatrix()[1:2, ]
  })
  output$rnaSM <- renderTable({
    rnaseqSampleMetadata()[1:2, ]
  })
  output$rnaSG <- renderText({
    input$rnaseqSampleGrouping
  })

  # Display the first 2 rows of the small RNAseq count matrix, sample metadata,
  # and the selected grouping variable
  output$smallRnaCM <- renderTable({
    smallRnaseqCountMatrix()[1:2, ]
  })
  output$smallRnaSM <- renderTable({
    smallRnaseqSampleMetadata()[1:3, ]
  })
  output$smallRnaSG <- renderText({
    input$smallRnaSampleGrouping
  })

  # Display the first 2 rows of the proteomics count matrix, sample metadata,
  # and the selected grouping variable
  output$proteomicsCM <- renderTable({
    proteomicsCountMatrix()[1:2, ]
  })
  output$proteomicsSM <- renderTable({
    proteomicsSampleMetadata()
  })
  output$proteomicsSG <- renderText({
    input$proteomicsSampleGrouping
  })

  # Display the first 2 rows of the ATACseq peaks for condition 1 and 2
  output$atacPeak1Test <- renderTable({
    atacPeak1()[1:2, ]
  })
  output$atacPeak2Test <- renderTable({
    atacPeak2()[1:2, ]
  })

  # Display the first 2 rows of the miRNA - target interactions
  output$miRNATargetTest <- renderTable({
    miRNATarget()[1:2, ]
  })

  # Display the first 2 rows of the TF - target interactions
  output$TFTargetTest <- renderTable({
    TFTarget()[1:2, ]
  })

  # Display the first 2 rows of the small RNA annotation
  output$smallRNAAnnotationTest <- renderTable({
    smallRNAAnnotation()[1:2, ]
  })

  # Display the first 2 rows of the protein name conversion
  output$proteinNameConversionTest <- renderTable({
    proteinNameConversion()[1:2, ]
  })

  # Display the genome assembly
  output$genAssemblyTest <- renderText({
    genAssembly()
  })
}



# Build the shiny app object
shinyApp(ui = ui, server = server)

# [END]
