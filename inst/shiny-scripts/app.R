# Purpose: Define the shiny app for the package
# Author: Jielin Yang
# Date: 2023-11-29
# Version: 1.0
# Bugs and Issues: None

# Load libraries
library(shiny)
library(shinyBS)  # for interactive elements
library(DT)  # for interactive tables

# Define some global variables
# Note that RNA, SMALLRNA, ATAC, and PROTEIN is already defined within the
# package
# The following variables are used to define external data types
# Note these definitions are only available within the app
MIRNA <- "miRNA_target"
TF <- "TF_target"
HSAPIENS <- "Human (hg38)"  # species + assembly to support multiple assemblies
MMUSCULUS <- "Mouse (mm10)"

# input data names
INPUT_LIST <- list(
  RNAseq = c("rnaseqCountMatrix", "rnaseqSampleMetadata",
             "rnaseqSampleGrouping"),
  smallRNAseq = c("smallRnaseqCountMatrix", "smallRnaseqSampleMetadata",
                  "smallRnaSampleGrouping"),
  proteomics = c("proteomicsCountMatrix", "proteomicsSampleMetadata",
                 "proteomicsSampleGrouping"),
  ATACseq = c("atacPeak1", "atacPeak2"),
  miRNA_target = "miRNATarget",
  TF_target = "TFTarget"
)

# Bioconductor annotation packages for different genome assemblies
BIOC_ANNOTATION <- list(
  "Human (hg38)" = list(
    organism = "org.Hs.eg.db",
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    bsgenome = "BSgenome.Hsapiens.UCSC.hg38"
  ),
  "Mouse (mm10)" = list(
    organism = "org.Mm.eg.db",
    txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    bsgenome = "BSgenome.Mmusculus.UCSC.mm10"
  )
)

# Options for handling data tables
# Input data
DT_OPTIONS_INPUT <- list(
      dom = "Brtip",
      scrollX = TRUE,
      scrollCollapse = TRUE
    )
DT_OPTIONS_DEOUTPUT <- list(
      pageLength = 10,
      lengthMenu = c(5, 10, 20, 50)
    )

#' Helper function to check for missing inputs
#'
#' @keywords internal
#'
#' @param input The input object from the shiny app for server code
#' @param requiredInputs A character vector of required input names
#'
#' @return A boolean value indicating whether there are missing inputs
#'
checkMissingInputs <- function(input, requiredInputs) {
  missing <- FALSE
  if (!all(requiredInputs %in% names(input))) {
    missing <- TRUE
  } else {
    for (i in requiredInputs) {
      if (is.null(input[[i]])) {
        missing <- TRUE
        break
      }
    }
  }
  return(missing)
}

#' Helper function to retrieve full DE results
#'
#' @keywords internal
#'
#' This function has a precondition that the DE analysis has already been
#' performed for the specified omic data type.
#'
#' @param objMOList A MOList object containing the DE results
#' @param omic A character string indicating the omic data type
#'
#' @return A data frame containing the full DE results
#'
getFullDEResults <- function(objMOList,
                             omic = c(RNA, SMALLRNA, ATAC, PROTEIN)) {
  omic <- match.arg(omic, c(RNA, SMALLRNA, ATAC, PROTEIN))
  deResult <- switch(omic,
                     "RNAseq" = objMOList$DERNAseq,
                     "smallRNAseq" = objMOList$DEsmallRNAseq,
                     "ATACseq" = objMOList$DEATAC,
                     "proteomics" = objMOList$DEproteomics)
  if (omic %in% COUNT_OMICS) {
    # Exploit the TOPTag class to retrieve results with order and pi-values
    deResult <- TOPTag(deResult,
                       logFCCutoff = 0,
                       pCutoff = 1,
                       topGenes = 1,
                       direction = "both") %>% exportDE()

  } else {
    # Use the PEAKTag as.data.frame generic
    deResult <- deResult %>% as.data.frame()
  }
  return(deResult)
}

#' Helper function to retrieve annotation databases based on genome assembly
#'
#' Depending on the genome assembly the user selected, this function will
#' retrieve the list of required annotation databases for the analysis. It will
#' check if the required annotation packages are installed, and install them if
#' not AT RUNTIME. Many of the annotation packages are large, so they are not
#' checked at package installation, but rather at runtime and are conditionally
#' installed upon user request. This function is only called immediately before
#' performing motif enrichment analysis to prevent unnecessary memory usage.
#'
#' @keywords internal
#'
#' @param genAssembly A character string indicating the genome assembly
#'
#' @return A list of annotation databases
getAnnotationDatabases <- function(genAssembly, session, input) {
  # Retrieve the list of annotation databases
  annotationDatabases <- BIOC_ANNOTATION[[genAssembly]]
  # Check if the annotation packages are installed
  stillMissing <- c()
  if (!requireNamespace(annotationDatabases$organism, quietly = TRUE)) {
    stillMissing <- c(stillMissing, annotationDatabases$organism)
  }
  if (!requireNamespace(annotationDatabases$txdb, quietly = TRUE)) {
    stillMissing <- c(stillMissing, annotationDatabases$txdb)
  }
  if (!requireNamespace(annotationDatabases$bsgenome, quietly = TRUE)) {
    stillMissing <- c(stillMissing, annotationDatabases$bsgenome)
  }
  # Display a modal to warn the user that the annotation packages are missing
  if (length(stillMissing) > 0) {
    showModal(
      modalDialog(
        title = paste0("Missing Annotation Packages for ", genAssembly),
        paste0(
          "The following annotation packages are missing: ",
          paste(stillMissing, collapse = ", "),
          ". Install them now?"
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("installAnno", "Install", class = "btn-primary")
        )
      )
    )
    # If the user clicks the install button, install the packages
    observeEvent(input$installAnno, {
      removeModal()
      showModal(modalDialog(
        title = "Installing missing annotation packages",
        paste0("Please wait while the packages are being installed... This ",
               "may take a while... ",
               "Please re-run differential analysis after the packages are ",
               "installed."),
      ))
      BiocManager::install(stillMissing, update = FALSE, ask = FALSE)
    })
    return(NULL)
  } else {
    # If the packages are already installed, retrieve the annotation databases
    # ChIPseeker requires loading org packages from the namespace
    library(annotationDatabases$organism, character.only = TRUE)
    # Retrieve the annotation databases with their names as variables without
    # loading the packages: solution retrieved from zkurtz at
    # https://stackoverflow.com/questions/48452555/
    # get-a-function-from-a-string-of-the-form-packagefunction?rq=3
    txdb <- getFromNamespace(annotationDatabases$txdb,
                              annotationDatabases$txdb)
    bsgenome <- getFromNamespace(annotationDatabases$bsgenome,
                                  annotationDatabases$bsgenome)
    return(list(organism = annotationDatabases$organism,
                txdb = txdb,
                bsgenome = bsgenome))
  }
}






# Although some dynamic UI code seems to be repetitive, it is necessary to
# define them separately to ensure that the UI elements are generated
# dynamically based on user inputs, and UI elements can only be generated
# within the UI/server code, so cannot refactor the code into a helper function


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

  # Title of the application
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
        ),

        # === Panel for section 4: differential analysis =======================
        tabPanel(
          "Section 4",
          tags$h4("Section 4: Exploring Differential Analysis"),
          br(),

          # Description of the section
          tags$p(tags$b("Description:"), "In this section, the app performs
          differential analysis of the multi-omics data to identify
          condition-specific transcriptomic alterations. Please provide
          settings to the program to define the behavior of the analysis."),
          br(),

          tags$h5("Please define the settings for the differential analysis:"),
          br(),
          tags$b("RNAseq differential expression:"),
          sliderInput(
            inputId = "rnaPadj",
            label = "Adjusted P-value threshold:",
            min = 0, max = 0.5, value = 0.05, step = 0.001
          ),
          sliderInput(
            inputId = "rnaLogFC",
            label = "Log2 fold change threshold:",
            min = 0, max = 10, value = 0, step = 0.1
          ),

          # Small RNAseq
          uiOutput("smallRnaCutoffs"),

          # Proteomics
          uiOutput("proteomicsCutoffs"),

          # ATACseq
          uiOutput("atacCutoffs"),

          # Batch effect correction variable selection
          br(),
          tags$h5("Please select the variable for batch effect correction:"),
          uiOutput("rnaBatchCorrection"),
          uiOutput("smallRnaBatchCorrection"),
          uiOutput("proteomicsBatchCorrection"),


          # === Starting differential analysis ===
          br(),
          tags$h5("Performing differential analysis:"),
          tags$p("By now, you should have uploaded all the required data
                and defined the settings for the differential analysis depending
                on the types of omic data selected."),
          # Method selection
          selectInput(
            inputId = "countDEMethod",
            label = tags$b("Select a method used for differential expression:"),
            choices = c(DESEQ2, EDGER),
            selected = DESEQ2
          ),
          # Button to start the analysis
          tags$p("Press the button below to
                 generate an initial set of differential analysis results."),
          bsButton(
            inputId = "runDifferentialAnalysis",
            label = "Run Differential Analysis",
            style = "primary"
          ),
          # Alert if not all required data are provided
          bsAlert(anchorId = "missingInputAlert"),
          # Alerts for error and warning messages during the analysis
          bsAlert(anchorId = "MOListErrorAlert"),
          bsAlert(anchorId = "diffOmicsErrorAlert"),
          bsAlert(anchorId = "annotateSmallRNAErrorAlert"),
          bsAlert(anchorId = "annotateProteinErrorAlert"),
          bsAlert(anchorId = "countPCAErrorAlert"),
          bsAlert(anchorId = "getAnnotationDatabasesErrorAlert"),
          bsAlert(anchorId = "diffMotifErrorAlert"),
          bsAlert(anchorId = "generateOutputErrorAlert"),

          # Updating results with different cutoffs without rerunning the
          # analysis
          br(),
          br(),
          tags$b("Update the cutoffs for the differential analysis:"),
          tags$p("If you would like to explore how different cutoffs affect
                the results, you can update the cutoffs without rerunning the
                analysis. The results will be updated automatically."),
        ),




        ),
    ),

    # Main panel for displaying outputs
    mainPanel(
      # Create a tabset panel for different output sections
      tabsetPanel(
        # All panels except for the first one are conditional panels, where
        # the panel is only shown if the corresponding data type is selected
        # and differential analysis is performed

        # === Panel 1: Input datasets ==========================================
        # Display the input data for users to double check their inputs

        tabPanel(
          "Data Input",
          tags$h4("Interactively Explore the Input Data"),
          br(),
          tags$p("The following sections display the raw input data provided. 
                 It is recommended to explore the data to ensure that the
                 correct data are provided. Raw data are displayed in tables
                 which you can sort and filter."),
          br(),
          tags$p(tags$b("Note:"), "Each data types are displayed in a separate
                panel. Depending on the data types selected, only titles of
                required data inputs are shown."),
          br(),
          # Make collapse panels since tables could be too large
          # Currently a bug in shinyBS prevents conditionally show/hide
          # collapse panels, with no workaround, so all panels are shown
          bsCollapse(
            bsCollapsePanel(
              title = "RNAseq Raw Data Input",
              tags$h5("Raw Count Matrix for RNAseq:"),
              DTOutput("rnaseqRawCountDT"),
              br(),
              tags$h5("Sample Metadata for RNAseq:"),
              DTOutput("rnaseqRawSampleMetadataDT")
            ),
            # Contents are conditionally shown for every data type except for
            # RNAseq
            bsCollapsePanel(
              title = "Small RNAseq Raw Data Input",
              uiOutput("smallRNACollapsePanel")
            ),
            bsCollapsePanel(
              title = "Proteomics Raw Data Input",
              uiOutput("proteomicsCollapsePanel")
            ),
            bsCollapsePanel(
              title = "ATACseq Raw Data Input",
              uiOutput("ATACseqCollapsePanel")
            ),
            bsCollapsePanel(
              title = "miRNA - Target Interactions",
              uiOutput("miRNATargetCollapsePanel")
            ),
            bsCollapsePanel(
              title = "TF - Target Interactions",
              uiOutput("TFTargetCollapsePanel")
            ),
            id = "collapseInput",
            open = "RNAseq Raw Data Input"
          ),
        ),

        # === Panel 2: RNAseq differential expression ==========================
        # Display the results of RNAseq differential expression analysis

        # display results only when DE analysis is performed
        tabPanel(
          "RNAseq Differential Expression",
          tags$h4("Differential Expression Analysis of RNAseq Data"),
          br(),

          # Plots
          tags$h5("Principle Component Analysis (PCA) of RNAseq Samples"),
          plotOutput(outputId = "rnaseqPCAPlot"),
          br(),
          tags$h5("Volcano Plot of RNAseq Data"),
          column(
            width = 12,
            align = "center",
            plotOutput(outputId = "rnaseqVolcanoPlot", width = "500px")
          ),
          br(),

          # Table of DE results
          tags$h5("Table of Differential Expression Results"),
          DTOutput(outputId = "rnaseqDETable"),
          tags$p(tags$b("Note: "),
                 "Selected rows will be highlighted in the volcano plot."),
          textOutput(outputId = "rnaseqGeneSelected"),
          br(),
          br(),
          tags$b("Download the analysis results:"),
          downloadButton(
            outputId = "downloadRnaseqDEResults",
            label = "Download RNAseq DE results (.csv)"
          ),
          downloadButton(
            outputId = "downloadRnaseqNormCounts",
            label = "Download RNAseq normalized counts (.csv)"
          ),
        ),

        # === Panel 3: Small RNAseq differential expression ====================
        tabPanel(
          "Small RNAseq Differential Expression",
          tags$h4("Differential Expression Analysis of Small RNAseq Data"),
          br(),

          # Plots
          tags$h5("Principle Component Analysis (PCA) of Small RNAseq Samples"),
          plotOutput(outputId = "smallRnaPCAPlot", height = "650px"),
          br(),
          tags$h5("Volcano Plot of Small RNAseq Data"),
          column(
            width = 6,
            align = "center",
            plotOutput(outputId = "smallRnaVolcanoPlotAnno", width = "500px"),
            tags$p("Color by small RNA type"),
          ),
          column(
            width = 6,
            align = "center",
            plotOutput(outputId = "smallRnaVolcanoPlotUpDown", width = "500px"),
            tags$p("Color by up/down regulation"),
          ),
          br(),

          # Table of DE results
          tags$h5("Table of Differential Expression Results"),
          DTOutput(outputId = "smallRnaDETable"),
          tags$p(tags$b("Note: "),
                 "Selected rows will be highlighted in the volcano plot."),
          textOutput(outputId = "smallRnaGeneSelected"),
          br(),
          br(),
          tags$b("Download the analysis results:"),
          downloadButton(
            outputId = "downloadSmallRnaDEResults",
            label = "Download Small RNAseq DE results (.csv)"
          ),
          downloadButton(
            outputId = "downloadSmallRnaNormCounts",
            label = "Download Small RNAseq normalized counts (.csv)"
          ),
        ),

        # === Panel 4: Proteomics differential expression ======================
        tabPanel(
          "Proteomics Differential Expression",
          tags$h4("Differential Expression Analysis of Proteomics Data"), 
          br(),

          # Plots
          tags$h5("Principle Component Analysis (PCA) of Proteomics Samples"),
          plotOutput(outputId = "proteomicsPCAPlot"),
          br(),
          tags$h5("Volcano Plot of Proteomics Data"),
          column(
            width = 12,
            align = "center",
            plotOutput(outputId = "proteomicsVolcanoPlot", width = "500px")
          ),
          br(),

          # Table of DE results
          tags$h5("Table of Differential Expression Results"),
          DTOutput(outputId = "proteomicsDETable"),
          tags$p(tags$b("Note: "),
                 "Selected rows will be highlighted in the volcano plot."),
          textOutput(outputId = "proteomicsGeneSelected"),
          br(),
          br(),
          tags$b("Download the analysis results:"),
          downloadButton(
            outputId = "downloadProteomicsDEResults",
            label = "Download Proteomics DE results (.csv)"
          ),
          downloadButton(
            outputId = "downloadProteomicsNormCounts",
            label = "Download Proteomics normalized counts (.csv)"
          ),
        ),

        # === Panel 5: ATACseq differential expression =========================
        tabPanel(
          "ATACseq Differential Motif Enrichment",
        ),










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

          # Test RNAseq cutoffs
          tags$p("RNAseq cutoffs:"),
          textOutput(outputId = "rnaPadjTest"),
          textOutput(outputId = "rnaLogFCTest"),

          # Test small RNAseq cutoffs
          tags$p("Small RNAseq cutoffs:"),
          textOutput(outputId = "smallRnaPadjTest"),
          textOutput(outputId = "smallRnaLogFCTest"),

          # Test proteomics cutoffs
          tags$p("Proteomics cutoffs:"),
          textOutput(outputId = "proteomicsPadjTest"),
          textOutput(outputId = "proteomicsLogFCTest"),

          # Test ATACseq cutoffs
          tags$p("ATACseq cutoffs:"),
          textOutput(outputId = "atacPadjTest"),
          textOutput(outputId = "atacUnadjPvalTest"),
          textOutput(outputId = "atacLogFCTest"),
          textOutput(outputId = "atacUseUnadjPTest"),

          # Test differential analysis method selection
          tags$p("Differential analysis method:"),
          textOutput(outputId = "deMethodTest"),



          # Show all existing inputs
          tags$p("All existing inputs:"),
          textOutput(outputId = "allInputs"),
          tags$p("Required inputs:"),
          textOutput(outputId = "requiredInputs"),

          # Print the MOList
          tags$p("MOList:"),
          verbatimTextOutput(outputId = "MOListTest"),

          # Show batch correction variable selection
          tags$p("Batch correction variable selection:"),
          textOutput(outputId = "rnaBatchVarTest"),
          textOutput(outputId = "smallRnaBatchVarTest"),
          textOutput(outputId = "proteomicsBatchVarTest"),

          # Run DE analysis?
          tags$p("Run DE analysis?"),
          textOutput(outputId = "runDETest"),

          # Printing databases
          tags$p("Annotation databases:"),
          verbatimTextOutput(outputId = "databasesTest"),



        ),
        selected = "Data Input",
        id = "mainOutputTabsetPanel"
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

  # Disable main output tabs until the analysis is performed
  hideTab(inputId = "mainOutputTabsetPanel", 
          target = "RNAseq Differential Expression")
  hideTab(inputId = "mainOutputTabsetPanel",
          target = "Small RNAseq Differential Expression")
  hideTab(inputId = "mainOutputTabsetPanel",
          target = "Proteomics Differential Expression")
  hideTab(inputId = "mainOutputTabsetPanel",
          target = "ATACseq Differential Motif Enrichment")

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
    if (useSmallRNAseq()) {
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
    } else {
      return(NULL)
    }
  })

  # Processing onsite download option
  output$onsiteDownloadSmallRnaseqCountMatrix <- downloadSmallRnaseqCountMatrix

  # Small RNAseq sample metadata (server)
  smallRnaseqSampleMetadata <- reactive({
    if (useSmallRNAseq()) {
      req(input$smallRnaseqSampleMetadata)
      # Validate file extension and read the file
      ext <- tools::file_ext(input$smallRnaseqSampleMetadata$datapath)
      validate(need(ext == "csv", "Please upload a .csv file"))
      smallRnaseqSampleMetadata <- read.csv(
        input$smallRnaseqSampleMetadata$datapath,
        header = TRUE, row.names = 1,
        stringsAsFactors = FALSE
      )
      return(smallRnaseqSampleMetadata)
    } else {
      return(NULL)
    }
  })

  # Processing onsite download option
  output$onsiteDownloadSmallRnaseqSampleMetadata <-
    downloadSmallRnaseqSampleMetadata

  # Retrieve the attribute (column) names of the sample metadata
  smallRnaseqSampleAttributes <- reactive({
    if (useSmallRNAseq()) {
      req(input$smallRnaseqSampleMetadata)
      colnames(smallRnaseqSampleMetadata())
    } else {
      return(NULL)
    }
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
    if (useProteomics()) {
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
    } else {
      return(NULL)
    }
  })

  # Processing onsite download option
  output$onsiteDownloadProteomicsCountMatrix <- downloadProteomicsCountMatrix

  # Proteomics sample metadata (server)
  proteomicsSampleMetadata <- reactive({
    if (useProteomics()) {
      req(input$proteomicsSampleMetadata)
      # Validate file extension and read the file
      ext <- tools::file_ext(input$proteomicsSampleMetadata$datapath)
      validate(need(ext == "csv", "Please upload a .csv file"))
      proteomicsSampleMetadata <- read.csv(
        input$proteomicsSampleMetadata$datapath,
        header = TRUE, row.names = 1,
        stringsAsFactors = FALSE
      )
      return(proteomicsSampleMetadata)
    } else {
      return(NULL)
    }
  })

  # Processing onsite download option
  output$onsiteDownloadProteomicsSampleMetadata <-
    downloadProteomicsSampleMetadata

  # Retrieve the attribute (column) names of the sample metadata
  proteomicsSampleAttributes <- reactive({
    if (useProteomics()) {
      req(input$proteomicsSampleMetadata)
      colnames(proteomicsSampleMetadata())
    } else {
      return(NULL)
    }
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
        tags$h5("Part ", partNumber, ": ATACseq Data Input"),
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
    if (useATACseq()) {
      req(input$atacPeak1)
      # Validate file extension and read the file
      ext <- tools::file_ext(input$atacPeak1$datapath)
      validate(need(ext == "bed", "Please upload a .bed file"))
      atacPeak1 <- read.table(input$atacPeak1$datapath, header = FALSE)
      atacPeak1 <- atacPeak1[, 1:3] %>% setNames(CHROMINFO)
      return(atacPeak1) # keep as data frame
    } else {
      return(NULL)
    }
  })

  # Processing onsite download option
  output$onsiteDownloadATACPeak1 <- downloadATACPeak1

  # ATACseq peaks for condition 2 (server)
  atacPeak2 <- reactive({
    if (useATACseq()) {
      req(input$atacPeak2)
      # Validate file extension and read the file
      ext <- tools::file_ext(input$atacPeak2$datapath)
      validate(need(ext == "bed", "Please upload a .bed file"))
      atacPeak2 <- read.table(input$atacPeak2$datapath, header = FALSE)
      atacPeak2 <- atacPeak2[, 1:3] %>% setNames(CHROMINFO)
      return(atacPeak2) # keep as data frame
    } else {
      return(NULL)
    }
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
    if (useSmallRNAseq()) {
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
    } else {
      return(NULL)
    }
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
    if (useProteomics()) {
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
    } else {
      return(NULL)
    }
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
          title = "Great! All small RNAs are annotated.",
          style = "success"
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
            "The protein-gene ID conversion infomration does not cover all ",
            "proteins in the ",
            "proteomics data. Uncovered proteins and their corresoinding ",
            "genes may be lost during TRN construction."
          )
        )
      } else {
        createAlert(session,
          anchorId = "proteinFullCoverageAlert",
          alertId = "proteinFullCoverageAlertID",
          title = "Great! All proteins are covered.",
          style = "success"
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
    if (usemiRNATarget()) {
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
    } else {
      return(NULL)
    }
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
    if (useTFTarget()) {
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
    } else {
      return(NULL)
    }
  })

  # Processing onsite download option
  output$onsiteDownloadTFTarget <- downloadTFTarget

  # Display the results of raw data input to the user in the Data Input
  # main panel
  # RNAseq
  output$rnaseqRawCountDT <- renderDT(
    rnaseqCountMatrix(),
    options = DT_OPTIONS_INPUT
  )
  output$rnaseqRawSampleMetadataDT <- renderDT(
    rnaseqSampleMetadata(),
    options = DT_OPTIONS_INPUT
  )

  # Small RNAseq
  output$smallRnaseqRawCountDT <- renderDT(
    smallRnaseqCountMatrix(),
    options = DT_OPTIONS_INPUT
  )
  output$smallRnaseqRawSampleMetadataDT <- renderDT(
    smallRnaseqSampleMetadata(),
    options = DT_OPTIONS_INPUT
  )
  output$smallRNACollapsePanel <- renderUI({
    if (useSmallRNAseq()) {
      tagList(
        tags$h5("Raw Count Matrix for Small RNAseq:"),
        DTOutput("smallRnaseqRawCountDT"),
        br(),
        tags$h5("Sample Metadata for RNAseq:"),
        DTOutput("smallRnaseqRawSampleMetadataDT"),
      )
    }
  })
  
  # Proteomics
  output$proteomicsRawCountDT <- renderDT(
    proteomicsCountMatrix(),
    options = DT_OPTIONS_INPUT
  )
  output$proteomicsRawSampleMetadataDT <- renderDT(
    proteomicsSampleMetadata(),
    options = DT_OPTIONS_INPUT
  )
  output$proteomicsCollapsePanel <- renderUI({
    if (useProteomics()) {
      tagList(
        tags$h5("Raw Count Matrix for Proteomics:"),
        DTOutput("proteomicsRawCountDT"),
        br(),
        tags$h5("Sample Metadata for Proteomics:"),
        DTOutput("proteomicsRawSampleMetadataDT"),
      )
    }
  })

  # ATACseq
  output$atacPeak1DT <- renderDT(
    atacPeak1(),
    options = DT_OPTIONS_INPUT
  )
  output$atacPeak2DT <- renderDT(
    atacPeak2(),
    options = DT_OPTIONS_INPUT
  )
  output$ATACseqCollapsePanel <- renderUI({
    if (useATACseq()) {
      tagList(
        tags$h5("Peaks for Condition 1:"),
        DTOutput("atacPeak1DT"),
        br(),
        tags$h5("Peaks for Condition 2:"),
        DTOutput("atacPeak2DT"),
      )
    }
  })

  # External miRNA - target interactions
  output$miRNATargetDT <- renderDT(
    miRNATarget(),
    options = DT_OPTIONS_INPUT
  )
  output$miRNATargetCollapsePanel <- renderUI({
    if (usemiRNATarget()) {
      tagList(
        tags$h5("miRNA - Target Interactions:"),
        DTOutput("miRNATargetDT"),
      )
    }
  })

  # External TF - target interactions
  output$TFTargetDT <- renderDT(
    TFTarget(),
    options = DT_OPTIONS_INPUT
  )
  output$TFTargetCollapsePanel <- renderUI({
    if (useTFTarget()) {
      tagList(
        tags$h5("TF - Target Interactions:"),
        DTOutput("TFTargetDT"),
      )
    }
  })

  # === Section 4: differential analysis =======================================

  # According to the data types selected, ask the user to define the cutoffs
  # for differential analysis

  # RNAseq cutoffs, defined in UI since mandatory
  rnaPadj <- reactive({
    req(input$rnaPadj)
    return(input$rnaPadj)
  })
  rnaLogFC <- reactive({
    req(input$rnaLogFC)
    return(input$rnaLogFC)
  })

  # small RNAseq cutoffs (dynamic UI)
  output$smallRnaCutoffs <- renderUI({
    if (useSmallRNAseq()) {
      tagList(
        # Small RNAseq p-value and logFC cutoffs
        br(),
        tags$b("Small RNAseq differential expression:"),
        sliderInput(
          inputId = "smallRnaPadj",
          label = "Adjusted P-value threshold",
          min = 0, max = 0.5, value = 0.05, step = 0.001
        ),
        sliderInput(
          inputId = "smallRnaLogFC",
          label = "Log2 fold change threshold:",
          min = 0, max = 10, value = 0, step = 0.1
        ),
      )
    }
  })
  # Server for data retrieval
  smallRnaPadj <- reactive({
    if (useSmallRNAseq()) {
      req(input$smallRnaPadj)
      return(input$smallRnaPadj)
    } else {
      return(NULL)
    }
  })
  smallRnaLogFC <- reactive({
    if (useSmallRNAseq()) {
      req(input$smallRnaLogFC)
      return(input$smallRnaLogFC)
    } else {
      return(NULL)
    }
  })

  # Proteomics cutoffs (dynamic UI)
  output$proteomicsCutoffs <- renderUI({
    if (useProteomics()) {
      tagList(
        # Proteomics p-value and logFC cutoffs
        br(),
        tags$b("Proteomics differential expression:"),
        sliderInput(
          inputId = "proteomicsPadj",
          label = "Adjusted P-value threshold",
          min = 0, max = 0.5, value = 0.05, step = 0.001
        ),
        sliderInput(
          inputId = "proteomicsLogFC",
          label = "Log2 fold change threshold:",
          min = 0, max = 10, value = 0, step = 0.1
        ),
      )
    }
  })
  # Server for data retrieval
  proteomicsPadj <- reactive({
    if (useProteomics()) {
      req(input$proteomicsPadj)
      return(input$proteomicsPadj)
    } else {
      return(NULL)
    }
  })
  proteomicsLogFC <- reactive({
    if (useProteomics()) {
      req(input$proteomicsLogFC)
      return(input$proteomicsLogFC)
    } else {
      return(NULL)
    }
  })

  # ATACseq cutoffs (dynamic UI)
  output$atacCutoffs <- renderUI({
    if (useATACseq()) {
      tagList(
        br(),
        # Checkbox to choose whether to use p-value cutoff
        tags$b("ATACseq differential motif enrichment:"),
        # Use conditionalPanel to hide/show one of the two options
        conditionalPanel(
          condition = "input.useATACunadjPval == false",
          sliderInput(
            inputId = "atacPadj",
            label = "Adjusted P-value threshold:",
            min = 0, max = 0.5, value = 0.05, step = 0.001
          )
        ),
        conditionalPanel(
          condition = "input.useATACunadjPval == true",
          sliderInput(
            inputId = "atacUnadjPval",
            label = "Unadjusted P-value threshold:",
            min = 0, max = 0.5, value = 0.05, step = 0.001
          )
        ),
        # A checkbox to choose whether to use unadjusted p-value
        checkboxInput(
          inputId = "useATACunadjPval",
          label = "Use unadjusted P-value instead?",
          value = FALSE
        ),
        div(style = "margin-top: -10px"),
        tags$p(id = "useATACunadjPvalNote", "Why using unadjusted P-value?
        (click here to see details)"),
        # add a pop-up to explain the use of unadjusted p-value
        bsPopover(
        id = "useATACunadjPvalNote",
        title = "Choosing P-value in differential motif enrichment analysis",
        content = paste0(
      "Differential motif enrichment analysis tests for transcription factor ",
      "or cis-regulatory motifs that are found to be differentially enriched ",
      "between two conditions. It requires that a motif must be enriched in ",
      "at least one condition, and at the same time there must be a ",
      "significant difference in enrichment between the two conditions. ",
      "Therefore, this analysis can be very stringent. It is up to the ",
      "users to decide whether to lessen the stringency by using unadjusted ",
      "p-value instead of adjusted p-value."), trigger = "click"),
        # and finally the fold enrichment cutoff
        sliderInput(
          inputId = "atacLogFC",
          label = "Log2 fold enrichment threshold:",
          min = 0, max = 10, value = 0, step = 0.1
        ),
      )
    }
  })

  # Server for data retrieval
  atacUseUnadjP <- reactive({
    return(input$useATACunadjPval)
  })
  atacPadj <- reactive({
    if (useATACseq()) {
      req(input$atacPadj)
      return(input$atacPadj)
    } else {
      return(NULL)
    }
  })
  atacUnadjPval <- reactive({
    if (useATACseq() && atacUseUnadjP()) {
      req(input$atacUnadjPval)
      return(input$atacUnadjPval)
    } else {
      return(NULL)  # return NULL if the checkbox is not checked
    }
  })
  atacLogFC <- reactive({
    if (useATACseq()) {
      req(input$atacLogFC)
      return(input$atacLogFC)
    } else {
      return(NULL)
    }
  })

  # Batch correction (dynamic UI)
  # Even though this input for RNAseq is always required, the presence of the
  # UI depends on whether user has uploaded RNAseq sample metadata
  output$rnaBatchCorrection <- renderUI({
    selectInput(
      inputId = "rnaBatchVar",
      label = "Batch variable for RNAseq:",
      choices = c("None", rnaseqSampleAttributes()),
      multiple = FALSE,
      selected = "NULL"
    )
  })
  output$smallRnaBatchCorrection <- renderUI({
    if (useSmallRNAseq()) {
      selectInput(
        inputId = "smallRnaBatchVar",
        label = "Batch variable for small RNAseq:",
        choices = c("None", smallRnaseqSampleAttributes()),
        multiple = FALSE,
        selected = "NULL"
      )
    }
  })
  output$proteomicsBatchCorrection <- renderUI({
    if (useProteomics()) {
      selectInput(
        inputId = "proteomicsBatchVar",
        label = "Batch variable for proteomics:",
        choices = c("None", proteomicsSampleAttributes()),
        multiple = FALSE,
        selected = "NULL"
      )
    }
  })

  # Server for batch variable retrieval
  rnaBatchVar <- reactive({
    return(input$rnaBatchVar)
  })                              # If user choose "None", returns NULL
  smallRnaBatchVar <- reactive({
    return(input$smallRnaBatchVar)
  })
  proteomicsBatchVar <- reactive({
    return(input$proteomicsBatchVar)
  })


  # === Perform differential analysis ==========================================

  # Retrieve DE method
  deMethod <- reactive({
    req(input$countDEMethod)
    return(input$countDEMethod)
  })

  rnaseqGeneIndex <- reactive({
    input$rnaseqDETable_rows_selected
  })

  # The analysis must follow a specific order of the pipeline
  # First need to initialize a reactive value to store the MOList object
  reactiveMOList <- reactiveValues(objMOList = NULL)
  # Logic for checking whether the analysis has been performed
  reactiveOverwrite <- reactiveValues(overwrite = FALSE, previous = FALSE)

  # Perform differential analysis on the button click
  observeEvent(input$runDifferentialAnalysis, {
    # First, close any alert messages if they are oped in prior runs
    closeAlert(session, "missingInputAlertID")
    closeAlert(session, "MOListErrorAlertID")
    closeAlert(session, "diffOmicsErrorAlertID")
    closeAlert(session, "diffOmicsWarningAlertID")
    closeAlert(session, "annotateSmallRNAErrorAlertID")
    closeAlert(session, "annotateProteinErrorAlertID")
    closeAlert(session, "countPCAErrorAlertID")
    closeAlert(session, "getAnnotationDatabasesErrorAlertID")
    closeAlert(session, "diffMotifErrorAlertID")

    # Initialize a variable to break the execution if any error occurs
    continueExecution <- TRUE

    # Depending on whether the analysis has been performed, create a modal
    # to warn users that re-running the analysis will overwrite the previous
    # results
    if (!is.null(reactiveMOList$objMOList)) {
      # DE analysis has been performed, open a modal
      reactiveOverwrite$previous <- TRUE
      showModal(
        modalDialog(
          title = "Re-run differential analysis?",
          "A previous differential analysis has been found. Re-running the ",
          "analysis will overwrite the previous results. Are you sure you ",
          "want to continue?",
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirmOverwrite", "Continue", class = "btn-primary")
          )
        )
      )
      # If the user confirms, set the reactive value to TRUE
      observeEvent(input$confirmOverwrite, {
        reactiveOverwrite$overwrite <- TRUE
        removeModal()
      })
    } else {
      # Do nothing
    }

    # Run the analysis only if this is a new analysis or the user confirms
    # overwriting the previous results
    if (!reactiveOverwrite$previous ||
        (reactiveOverwrite$previous && reactiveOverwrite$overwrite)) {
      # Perform the analysis
      # Create a Progress object
      progress <- shiny::Progress$new(max = 10)
      # Make sure it closes when exiting the function
      on.exit(progress$close())

      # Step 1: Construct the MOList object
      progress$set(message = "Constructing the MOList object...", value = 1)
      # If any input is missing, stop the analysis
      requiredInputs <- INPUT_LIST[omicsDataTypes()] %>% unlist()
      if (useSmallRNAseq()) {
        requiredInputs <- c(requiredInputs, "smallRNAAnnotation")
      }
      if (useProteomics()) {
        requiredInputs <- c(requiredInputs, "proteinNameConversion")
      }
      # Check all required inputs exist
      if (checkMissingInputs(input, requiredInputs)) {
        continueExecution <<- FALSE
        # Render a warning message
        createAlert(session,
          anchorId = "missingInputAlert",
          alertId = "missingInputAlertID",
          title = "Missing Input",
          content = paste0(
            "Please make sure all required inputs are provided in the ",
            "previous sections."
          ),
          style = "danger"
        )
      } else {
        # Continue
      }

      # Proceed to construct the MOList object
      # Internal error checking of the functions are used, but they will cause
      # the app to crash. Must be wrapped in tryCatch block to handle
      # exceptions
      objMOList <- NULL   # must initialize outside of tryCatch

      if (continueExecution){
        tryCatch({
          # Construct the MOList object
          objMOList <- MOList(
            RNAseq = rnaseqCountMatrix(),
            RNAGroupBy = rnaseqSampleMetadata()[[input$rnaseqSampleGrouping]],
            smallRNAseq = smallRnaseqCountMatrix(),
            smallRNAGroupBy =
              smallRnaseqSampleMetadata()[[input$smallRnaSampleGrouping]],
            proteomics = proteomicsCountMatrix(),
            proteomicsGroupBy =
              proteomicsSampleMetadata()[[input$proteomicsSampleGrouping]],
            ATACpeak1 = atacPeak1(),
            ATACpeak2 = atacPeak2()
          )
        }, error = function(e) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "MOListErrorAlert",
            alertId = "MOListErrorAlertID",
            title = "Error in handling omic input data:",
            content = paste0(
              "<strong> Please make sure the input data are in the correct ",
              "format. </strong><br>",
              "<strong> Ecountered the following error: </strong><br>",
              conditionMessage(e)
            ),
            style = "danger"
          )
        })
      }

      # Step 2: Perform differential analysis
      omicsData <- setdiff(omicsDataTypes(), c(MIRNA, TF))
      progress$set(message = "Performing differential analysis...",
                    detail = paste0("on ", paste(omicsData, collapse = ", ")),
                    value = 4)
      if (continueExecution){
        tryCatch({
          objMOList <- diffOmics(objMOList,
            rnaseqBatch = rnaseqSampleMetadata()[[rnaBatchVar()]],
            smallRnaBatch =
              smallRnaseqSampleMetadata()[[smallRnaBatchVar()]],
            proteinBatch =
              proteomicsSampleMetadata()[[proteomicsBatchVar()]],
            program = deMethod()
          )
        }, error = function(e) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "diffOmicsErrorAlert",
            alertId = "diffOmicsErrorAlertID",
            title = "Error in performing differential analysis:",
            content = paste0(
              "<strong> Please make sure the input data are valid for the ",
              "selected differential analysis method. </strong><br>",
              "<strong> Ecountered the following error: </strong><br>",
              conditionMessage(e)
            ),
            style = "danger"
          )
        })
      }

      # Step 3: Adding annotations if required
      progress$set(message = "Checking for annotations...", value = 5)
      if (continueExecution){
        tryCatch({
          if (useSmallRNAseq()) {
            objMOList <- annotateSmallRNA(objMOList,
                                          anno = smallRNAAnnotation())
          } else {
            # Continue
          }
        }, error = function(e) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "annotateSmallRNAErrorAlert",
            alertId = "annotateSmallRNAErrorAlertID",
            title = "Error in annotating small RNAseq data:",
            content = paste0(
              "<strong> Please make sure the annotation is provided ",
              "in the correct format </strong><br>",
              "<strong> Ecountered the following error: </strong><br>",
              conditionMessage(e)
            ),
            style = "danger"
          )
        })
      }
      if (continueExecution){
        tryCatch({
          if (useProteomics()) {
            objMOList <- setGene2Protein(objMOList,
                                         proteinNameConversion())
          } else {
            # Continue
          }
        }, error = function(e) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "annotateProteinErrorAlert",
            alertId = "annotateProteinErrorAlertID",
            title = "Error in annotating proteomics data:",
            content = paste0(
              "<strong> Please make sure the annotation is provided ",
              "in the correct format </strong><br>",
              "<strong> Ecountered the following error: </strong><br>",
              conditionMessage(e)
            ),
            style = "danger"
          )
        })
      }

      # Step 4: Principal component analysis
      # RNAseq
      progress$set(message = "Performing principal component analysis...",
                    detail = paste0("on ",
                            paste(setdiff(omicsData, ATAC), collapse = ", ")),
                    value = 6)
      pltRNAPCA <- NULL   # Again, must initialize outside of tryCatch
      pltProteomicsPCA <- NULL
      pltSmallRNAPCA <- NULL
      if (continueExecution){
        tryCatch({
          pltRNAPCA <- countPCA(matrix = getRawData(objMOList, RNA),
                                groupBy = objMOList@RNAseqSamples$groupBy,
                                batch = rnaseqSampleMetadata()[[rnaBatchVar()]])
          # Proteomics
          if (useProteomics()) {
          pltProteomicsPCA <- countPCA(matrix = getRawData(objMOList, PROTEIN),
                    groupBy = objMOList@proteomicsSamples$groupBy,
                    batch = proteomicsSampleMetadata()[[proteomicsBatchVar()]])
          } else {
            # Continue
          }
          # Small RNAseq
          if (useSmallRNAseq()) {
            pltSmallRNAPCA <- plotSmallRNAPCAs(objMOList,
                      batch = smallRnaseqSampleMetadata()[[smallRnaBatchVar()]])
          } else {
            # Continue
          }
        }, error = function(e) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "countPCAErrorAlert",
            alertId = "countPCAErrorAlertID",
            title = "Error in performing principal component analysis:",
            content = paste0(
              "<strong> Ecountered the following error: </strong><br>",
              conditionMessage(e)
            ),
            style = "danger"
          )
        })
      }

      # Step 5: Differential motif enrichment analysis
      progress$set(message = "Performing differential motif enrichment
                              analysis...",
                   detail = "Assigning features to peaks... Searching for
                             motifs... This may take a while...",
                   value = 8)
      # Retrieve the annotation databases
      annoDB <- NULL
      if (useATACseq()) {
        annoDB <- getAnnotationDatabases(genAssembly(), session, input)
        if (is.null(annoDB)) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "getAnnotationDatabasesErrorAlert",
            alertId = "getAnnotationDatabasesErrorAlertID",
            title = "Error in retrieving annotation databases:",
            content = paste0(
              "<strong> No annotation databases found. </strong><br>"
            ),
            style = "danger"
          )
        } else {
          # Continue
        }
      } else {
        # Continue
      }
      # Perform the Analysis
      if (continueExecution && useATACseq() && !is.null(annoDB)) {
        tryCatch({
          objMOList <- annotateATACPeaksMotif(
            objMOList,
            TxDb = annoDB$txdb,
            annoDb = annoDB$organism,
            bsgenome = annoDB$bsgenome,
            pwmL = jasparVertebratePWM
          )
        }, error = function(e) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "diffMotifErrorAlert",
            alertId = "diffMotifErrorAlertID",
            title = "Error in performing differential motif enrichment
                    analysis:",
            content = paste0(
              "<strong> Ecountered the following error: </strong><br>",
              conditionMessage(e)
            ),
            style = "danger"
          )
        })
      } else {
        # Continue
      }

      # Step 6: Generate the output
      progress$set(message = "Generating final outputs...",
                   detail = paste0("for ",
                                   paste(omicsDataTypes(), collapse = ", ")),
                   value = 9)
      if (continueExecution){
        tryCatch({
          # Render plots & tables to the main panel
          # RNAseq DE table
          # Options used for table obtained from DT package documentation:
          # https://rstudio.github.io/DT/options.html
          rnaDF <- getFullDEResults(objMOList, omic = RNA)
          # Show only 3 digits with scientific notation
          rnaDF[, -1] <- format(rnaDF[, -1], digits = 4, scientific = TRUE)
          output$rnaseqDETable <- renderDT(
            rnaDF,
            filter = "top",
            options = DT_OPTIONS_DEOUTPUT
          )
          # RNAseq volcano plot
          # Plot it
          output$rnaseqVolcanoPlot <- renderPlot({
            plotVolcano(objMOList, omic = RNA,
                        adjP = rnaPadj(), log2FC = rnaLogFC(), 
                        highlight = rownames(rnaDF)[rnaseqGeneIndex()])
          })
          # RNAseq PCA plot
          output$rnaseqPCAPlot <- renderPlot({
            pltRNAPCA
          })
          # Selected Genes
          output$rnaseqGeneSelected <- renderText({
            if (length(rnaseqGeneIndex()) > 0) {
              paste0("Gene selected: ", 
                paste(rownames(rnaDF)[rnaseqGeneIndex()], collapse = ", "))
            } else {
              # Do nothing
            }
          })

          # Small RNAseq 
          if (useSmallRNAseq()) {
            # Small RNAseq DE table
            smallRnaDF <- getFullDEResults(objMOList, omic = SMALLRNA)
            # Show only 3 digits with scientific notation
            smallRnaDF[, -1] <- format(smallRnaDF[, -1], digits = 4,
                                       scientific = TRUE)
            output$smallRnaDETable <- renderDT(
              smallRnaDF,
              filter = "top",
              options = DT_OPTIONS_DEOUTPUT
            )
            # Small RNAseq volcano plot
            output$smallRnaVolcanoPlotAnno <- renderPlot({
              plotVolcanoSmallRNA(objMOList, adjP = smallRnaPadj(),
                                  log2FC = smallRnaLogFC(),
                                  highlight = rownames(smallRnaDF)[
                                    input$smallRnaDETable_rows_selected])
            })
            output$smallRnaVolcanoPlotUpDown <- renderPlot({
              plotVolcano(objMOList, omic = SMALLRNA,
                          adjP = smallRnaPadj(), log2FC = smallRnaLogFC(),
                          highlight = rownames(smallRnaDF)[
                            input$smallRnaDETable_rows_selected])
            })
            # Small RNAseq PCA plots
            output$smallRnaPCAPlot <- renderPlot({
              pltSmallRNAPCA$miRNA + pltSmallRNAPCA$piRNA +
                pltSmallRNAPCA$circRNA + pltSmallRNAPCA$snoRNA +
                pltSmallRNAPCA$snRNA + pltSmallRNAPCA$tRNA
            })
            # Selected genes
            output$smallRnaGeneSelected <- renderText({
              if (length(input$smallRnaDETable_rows_selected) > 0) {
                paste0("Gene selected: ", 
                  paste(rownames(smallRnaDF)[
                    input$smallRnaDETable_rows_selected], collapse = ", "))
              } else {
                # Do nothing
              }
            })
          } else {
            # Continue
          }

          # Proteomics
          if (useProteomics()) {
            # Proteomics DE table
            proteomicsDF <- getFullDEResults(objMOList, omic = PROTEIN)
            # Show only 3 digits with scientific notation
            proteomicsDF[, -1] <- format(proteomicsDF[, -1], digits = 4,
                                         scientific = TRUE)
            output$proteomicsDETable <- renderDT(
              proteomicsDF,
              filter = "top",
              options = DT_OPTIONS_DEOUTPUT
            )
            # Proteomics volcano plot
            output$proteomicsVolcanoPlot <- renderPlot({
              plotVolcano(objMOList, omic = PROTEIN,
                          adjP = proteomicsPadj(), log2FC = proteomicsLogFC(),
                          highlight = rownames(proteomicsDF)[
                            input$proteomicsDETable_rows_selected])
            })
            # Proteomics PCA plot
            output$proteomicsPCAPlot <- renderPlot({
              pltProteomicsPCA
            })
            # Selected genes
            output$proteomicsGeneSelected <- renderText({
              if (length(input$proteomicsDETable_rows_selected) > 0) {
                paste0("Gene selected: ", 
                  paste(rownames(proteomicsDF)[
                    input$proteomicsDETable_rows_selected], collapse = ", "))
              } else {
                # Do nothing
              }
            })
          } else {
            # Continue
          }

          # ATACseq
          if (useATACseq()) {
            # ATACseq motif enrichment table
            atacDF <- getFullDEResults(objMOList, omic = ATAC)
            output$atacMotifTable <- renderDT(
              atacDF,
              filter = "top",
              options = c(DT_OPTIONS_DEOUTPUT, scrollX = TRUE)
            )
            # ATACseq annotation plot
            output$atacAnnotationPlot <- renderPlot({
              plotATACAnno(objMOList)
            })
            # ATACseq coverage plot
            output$atacCoveragePlot <- renderPlot({
              plotATACCoverage(objMOList)
            })
            # ATACseq motif enrichment heatmap

            output$atacMotifHeatmap <- renderPlot({
              plotATACMotifHeatmap(objMOList,
                                  pValueAdj = atacPadj(),
                                  pValue = atacUnadjPval(),
                                  log2FEnrich = atacLogFC())
            })
          } else {
            # Continue
          }

          # Generate the output
          reactiveMOList$objMOList <- objMOList

          # The tabPanels for each of the data types are generated dynamically
          # Also render download buttons for each of the data types
          showTab(inputId = "mainOutputTabsetPanel", 
                  target = "RNAseq Differential Expression",
                  select = TRUE)
          output$downloadRnaseqDEResults <- downloadHandler(
            filename = function() {
              paste("RNAseq_DE_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
              write.csv(getFullDEResults(objMOList, omic = RNA), file,
                        row.names = TRUE)
            }
          )
          output$downloadRnaseqNormCounts <- downloadHandler(
            filename = function() {
              paste("RNAseq_normalized_counts_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
              write.csv(exportNormalizedCounts(objMOList$DERNAseq), file,
                        row.names = TRUE)
            }
          )

          if (useSmallRNAseq()) {
            showTab(inputId = "mainOutputTabsetPanel", 
                    target = "Small RNAseq Differential Expression")
            output$downloadSmallRnaDEResults <- downloadHandler(
              filename = function() {
                paste("Small_RNAseq_DE_results_", Sys.Date(), ".csv", sep = "")
              },
              content = function(file) {
                write.csv(getFullDEResults(objMOList, omic = SMALLRNA), file,
                          row.names = TRUE)
              }
            )
            output$downloadSmallRnaNormCounts <- downloadHandler(
              filename = function() {
                paste("Small_RNAseq_normalized_counts_", Sys.Date(), ".csv",
                      sep = "")
              },
              content = function(file) {
                write.csv(exportNormalizedCounts(objMOList$DEsmallRNAseq),
                          file, row.names = TRUE)
              }
            )
          } else {
            # Continue
          }
          if (useProteomics()) {
            showTab(inputId = "mainOutputTabsetPanel", 
                    target = "Proteomics Differential Expression")
            output$downloadProteomicsDEResults <- downloadHandler(
              filename = function() {
                paste("Proteomics_DE_results_", Sys.Date(), ".csv", sep = "")
              },
              content = function(file) {
                write.csv(getFullDEResults(objMOList, omic = PROTEIN), file,
                          row.names = TRUE)
              }
            )
            output$downloadProteomicsNormCounts <- downloadHandler(
              filename = function() {
                paste("Proteomics_normalized_counts_", Sys.Date(), ".csv",
                      sep = "")
              },
              content = function(file) {
                write.csv(exportNormalizedCounts(objMOList$DEproteomics),
                          file, row.names = TRUE)
              }
            )
          } else {
            # Continue
          }
          if (useATACseq()) {
            showTab(inputId = "mainOutputTabsetPanel", 
                    target = "ATACseq Differential Motif Enrichment")
            output$downloadAtacMotifResults <- downloadHandler(
              filename = function() {
                paste("ATACseq_annotated _peaks_", Sys.Date(), ".csv",
                      sep = "")
              },
              content = function(file) {
                write.csv(getFullDEResults(objMOList, omic = ATAC), file,
                          row.names = TRUE)
              }
            )
          } else {
            # Continue
          }
        }, error = function(e) {
          continueExecution <<- FALSE
          createAlert(session,
            anchorId = "generateOutputErrorAlert",
            alertId = "generateOutputErrorAlertID",
            title = "Error in generating final outputs:",
            content = paste0(
              "<strong> Ecountered the following error: </strong><br>",
              conditionMessage(e)
            ),
            style = "danger"
          )
        })
      } else {
        # Continue
      }
      progress$set(message = "Done!", detail = NULL, value = 10)
    } 
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

  # Display the RNAseq cutoffs
  output$rnaPadjTest <- renderText({
    rnaPadj()
  })
  output$rnaLogFCTest <- renderText({
    rnaLogFC()
  })

  # Display the small RNAseq cutoffs
  output$smallRnaPadjTest <- renderText({
    smallRnaPadj()
  })
  output$smallRnaLogFCTest <- renderText({
    smallRnaLogFC()
  })

  # Display the proteomics cutoffs
  output$proteomicsPadjTest <- renderText({
    proteomicsPadj()
  })
  output$proteomicsLogFCTest <- renderText({
    proteomicsLogFC()
  })

  # Display the ATACseq cutoffs
  output$atacPadjTest <- renderText({
    atacPadj()
  })
  output$atacUnadjPvalTest <- renderText({
    atacUnadjPval()
  })
  output$atacLogFCTest <- renderText({
    atacLogFC()
  })
  output$atacUseUnadjPTest <- renderText({
    atacUseUnadjP()
  })

  # Differential analysis method
  output$deMethodTest <- renderText({
    deMethod()
  })

  # Showing all existing inputs
  output$allInputs <- renderText({
    paste(names(input), collapse = ", ")
  })

  # Show all required inputs
  reqInputs <- reactive({
    requiredInputs <- INPUT_LIST[omicsDataTypes()] %>% unlist()
    if (useSmallRNAseq()) {
      requiredInputs <- c(requiredInputs, "smallRNAAnnotation")
    }
    if (useProteomics()) {
      requiredInputs <- c(requiredInputs, "proteinNameConversion")
    }
    return(requiredInputs)
  })
  output$requiredInputs <- renderText({
    paste(reqInputs(), collapse = ", ")
  })

  # Show batch variables
  output$rnaBatchVarTest <- renderText({
    rnaBatchVar()
  })
  output$smallRnaBatchVarTest <- renderText({
    smallRnaBatchVar()
  })
  output$proteomicsBatchVarTest <- renderText({
    proteomicsBatchVar()
  })

  # Test if DE analysis should be run
  output$runDETest <- renderText({
    run <- !reactiveOverwrite$previous ||
        (reactiveOverwrite$previous && reactiveOverwrite$overwrite)
    return(as.character(run))
  })





}



# Build the shiny app object
shinyApp(ui = ui, server = server)

# [END]
