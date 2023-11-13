
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IntegraTRN

------------------------------------------------------------------------

<!-- badges: start -->

[![License: GPL (\>=
3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://choosealicense.com/licenses/gpl-3.0/)
![GitHub issues](https://img.shields.io/github/issues/j-y26/IntegraTRN)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Description

The R package `IntegraTRN` integrates transcriptomic, small RNAomic,
proteomic, and epigenomic data to construct a transcriptional regulatory
network (TRN) underlying the gene expression changes in a developmental
or disease (or any continuous/binary) biological context. In particular,
the TRN is a network of interacting factors, in which one element can
exert a regulatory effect (activating or inhibitory) on the expression
of one or more elements. Given the complex nature of transcriptional
regulation, the package reveals the key players of such regulation,
including transcriptional factors (TFs) and small RNAs. Since
gene/protein expression changes are usually the most direct molecular
cause of observed phenotypic alterations, the core TRN deciphers changes
between two or more conditions and infers the upstream regulatory
factors that mediate such changes. The analysis pipeline provided by
this package is primarily composed of a two-step process: (1)
elucidating the transcriptional and regulatory changes that take place
during a biological process, such as development or disease; and (2)
performing correlational analysis with rigorous filtering based on
biological data to identify the regulatory interactions that are
responsible for the observed changes.

The package is developed under the following environment:

- R version 4.3.1 (2023-06-16 ucrt)
- Platform: x86_64-w64-mingw32/x64 (64-bit)
- Running under: Windows 11 x64 (build 22621)

## Installation

To install the latest developmental version of the package:

``` r
require("devtools")
devtools::install_github("j-y26/IntegraTRN", build_vignettes = TRUE)
library("IntegraTRN")
```

To run the shinyApp: `Under construction`

## Overview

``` r
ls("package:IntegraTRN")
data(package = "IntegraTRN")
browseVignettes("IntegraTRN")
```

To get started, the functionality of the package is primarily divided
into two parts: (1) exploring differential expression or accessibility
of genes and proteins, and (2) constructing a small RNA - transcription
factor - gene regulatory network. The following is an overview of the
functions in the package, in the order of use subjected to the types of
omcis data available. `IntegraTRN` provides functions:

1.  ***MOList*** for generating a `MOList` object that contains the
    omics data and sample grouping information
2.  ***diffOmics*** for performing differential analysis on the omics
    data
3.  ***annotateSmallRNA*** for annotating small RNA transcripts
4.  ***plotVolcanoRNA*** for visualizing RNAseq differential expression
5.  ***plotVolcanoSmallRNA*** for visualizing small RNAseq differential
    expression
6.  ***plotSmallRNAPCAs*** for visualizing small RNAseq principal
    component analysis for each type of small RNA
7.  ***annotateATACPeaksMotif*** for annotating ATACseq peaks with motif
    enrichment analysis
8.  ***plotATACAnno*** for visualizing the annotation of ATACseq peaks
9.  ***plotATACCoverage*** for visualizing the coverage of ATACseq peaks
10. ***plotATACMotifHeatmap*** for visualizing the motif enrichment
    analysis as a heatmap comparing differentially enriched motifs
    between the two testing conditions
11. ***matchSamplesRNAsmallRNA*** for matching the samples between
    RNAseq and small RNAseq data
12. ***exportMatchResult*** for exporting the matching results to a data
    frame
13. ***loadExtInteractions*** for loading external interaction data for
    small RNA - gene and TF - gene interactions
14. ***setOmicCutoffs*** for setting the cutoffs for differential
    expression and accessibility used to filter the key elements in the
    TRN
15. ***constructTRN*** for constructing the TRN
16. ***plotNetwork*** for visualizing the TRN

The package also provides several datasets: - An RNAseq count matrix

- A small RNAseq count matrix coveting miRNA, tRNA, piRNA, snoRNA,
  snRNA, and circRNA

- A proteomics count matrix

- Sample information for all above 3 omics data

- Two ATACseq peak files as raw data located in the `extdata` folder

- An example miRNA-gene interaction dataset

- An example TF-gene interaction dataset

- An example MOList object containing all types of omics data, but with
  a very light weight (100 genes only)

Please refer to the package vignette for more details on these datasets
illustrates the analysis pipeline of the package.

An overview of the analysis pipeline is shown below:

<figure>
<img src="./inst/extdata/Schematics.jpg"
alt="Overview of the package analysis pipeline" />
<figcaption aria-hidden="true">Overview of the package analysis
pipeline</figcaption>
</figure>

## Contributions

The author of the package is Jielin Yang. The author defined all data
structures used in this package, including the S4 classes `MOList`,
`DETag`, `TOPTag`, `PEAKTag`, and `TRNet`. The author wrote the `MOList`
function to construct the key data structure and performs validations on
the input omics data. The package `GenomicTools.fileHandler` is used to
read bed files. The function `diffOmics` performs differential
expression on RNAseq, small RNAseq, and proteomics data using a negative
binomial model, which internally normalizes and performs differential
analysis using the `DESeq2` or `edgeR` package. ATACseq peaks are
handled as genomic coordinates using the `GenomicRanges` package. The
author wrote the `annotateSmallRNA` function to annotate small RNA
transcripts. The three plotting functions, `plotVolcanoRNA`,
`plotVolcanoSmallRNA`, and `plotSmallRNAPCAs`, are supported by the
`ggplot2` package, with PCA analysis supported by `DESeq2` on normalized
expression. The `ChIPseeker` package is used to annotate ATACseq peaks
with motif enrichment analysis and performs plotting on ATACseq peaks.
The `matchSamplesRNAsmallRNA` function performs optimal matching based
on mahalanobis distance between the sample information of the RNAseq and
small RNAseq data. The calculation of the mahalanobis distance and
selection of optimal pairs are supported by the `MatchIt` package. The
author wrote the utility functions `exportMatchResult` and
`loadExtInteractions` to export the matching results and load external
interaction data. The author also wrote the `setOmicCutoffs` function as
an easy way for the users to decide on the inclusion criteria for the
key elements in the TRN. The author wrote the `constructTRN` function,
with a logic defined to integrate the different omics data depending on
their availability, as well as whether predicted inference is used. The
inference of predicted small RNA - gene interactions is supported
partially by the author’s discretion to generate a single coherent
normalized expression matrix for both RNAseq and small RNAseq data that
allows co-expression estimation. The inference of small RNA - gene
interactions is performed by the `GENIE3` package, which internally uses
a three-based algorithm to infer the interactions. The `igraph` package
is used to visualize the TRN. Most data frame processing used internally
in the functions is supported by the `dplyr` package. Generative AI tool
was used to generate some unit test example data based on the author’s
description. Generative AI results were incorporated into the tests at
the author’s discretion.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-adar2023integrate" class="csl-entry">

Aharon-Yariv, Adar, Yaxu Wang, Abdalla Ahmed, and Paul Delgado-Olguin.
2023. “Integrated Small RNA, mRNA and Protein Omics Reveal a miRNA
Network Orchestrating Metabolic Maturation of the Developing Human
Heart.” *BMC Genomics*.

</div>

<div id="ref-ahmed2021maternal" class="csl-entry">

Ahmed, Abdalla, Minggao Liang, Lijun Chi, Yu-Qing Zhou, John G Sled,
Michael D Wilson, and Paul Delgado-Olguı́n. 2021. “Maternal Obesity
Persistently Alters Cardiac Progenitor Gene Expression and Programs
Adult-Onset Heart Disease Susceptibility.” *Molecular Metabolism* 43:
101116.

</div>

<div id="ref-networkd3" class="csl-entry">

Allaire, JJ, Peter Ellis, Christopher Gandrud, Kevin Kuo, BW Lewis,
Jonathan Owen, Kenton Russell, et al. 2017. *networkD3: D3 JavaScript
Network Graphs from r*. <https://CRAN.R-project.org/package=networkD3>.

</div>

<div id="ref-chang2020mirnet" class="csl-entry">

Chang, Le, Guangyan Zhou, Othman Soufan, and Jianguo Xia. 2020. “miRNet
2.0: Network-Based Visual Analytics for miRNA Functional Analysis and
Systems Biology.” *Nucleic Acids Research* 48 (W1): W244–51.

</div>

<div id="ref-csardi2006igraph" class="csl-entry">

Csardi, Gabor, Tamas Nepusz, et al. 2006. “The Igraph Software Package
for Complex Network Research.” *InterJournal, Complex Systems* 1695 (5):
1–9.

</div>

<div id="ref-genomictools" class="csl-entry">

Fischer, Daniel. 2020. *GenomicTools.fileHandler: File Handlers for
Genomic Data Analysis*.
<https://CRAN.R-project.org/package=GenomicTools.fileHandler>.

</div>

<div id="ref-huynh2010inferring" class="csl-entry">

Huynh-Thu, Vân Anh, Alexandre Irrthum, Louis Wehenkel, and Pierre
Geurts. 2010. “Inferring Regulatory Networks from Expression Data Using
Tree-Based Methods.” *PloS One* 5 (9): e12776.

</div>

<div id="ref-lawrence2013software" class="csl-entry">

Lawrence, Michael, Wolfgang Huber, Hervé Pagès, Patrick Aboyoun, Marc
Carlson, Robert Gentleman, Martin T Morgan, and Vincent J Carey. 2013.
“Software for Computing and Annotating Genomic Ranges.” *PLoS
Computational Biology* 9 (8): e1003118.

</div>

<div id="ref-li2020compsra" class="csl-entry">

Li, Jiang, Alvin T Kho, Robert P Chase, Lorena Pantano, Leanna Farnam,
Sami S Amr, and Kelan G Tantisira. 2020. “COMPSRA: A COMprehensive
Platform for Small RNA-Seq Data Analysis.” *Scientific Reports* 10 (1):
4552.

</div>

<div id="ref-love2014moderated" class="csl-entry">

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 1–21.

</div>

<div id="ref-machlab2022monalisa" class="csl-entry">

Machlab, Dania, Lukas Burger, Charlotte Soneson, Filippo M Rijli, Dirk
Schübeler, and Michael B Stadler. 2022. “monaLisa: An r/Bioconductor
Package for Identifying Regulatory Motifs.” *Bioinformatics* 38 (9):
2624–25.

</div>

<div id="ref-robinson2010edger" class="csl-entry">

Robinson, Mark D, Davis J McCarthy, and Gordon K Smyth. 2010. “edgeR: A
Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data.” *Bioinformatics* 26 (1): 139–40.

</div>

<div id="ref-stuart2011matchit" class="csl-entry">

Stuart, Elizabeth A, Gary King, Kosuke Imai, and Daniel Ho. 2011.
“MatchIt: Nonparametric Preprocessing for Parametric Causal Inference.”
*Journal of Statistical Software*.

</div>

<div id="ref-villanueva2019ggplot2" class="csl-entry">

Villanueva, Randle Aaron M, and Zhuo Job Chen. 2019. “Ggplot2: Elegant
Graphics for Data Analysis.” Taylor & Francis.

</div>

<div id="ref-dplyr" class="csl-entry">

Wickham, Hadley, Romain François, Lionel Henry, Kirill Müller, and Davis
Vaughan. 2023. *Dplyr: A Grammar of Data Manipulation*.
<https://CRAN.R-project.org/package=dplyr>.

</div>

<div id="ref-yu2015chipseeker" class="csl-entry">

Yu, Guangchuang, Li-Gen Wang, and Qing-Yu He. 2015. “ChIPseeker: An
r/Bioconductor Package for ChIP Peak Annotation, Comparison and
Visualization.” *Bioinformatics* 31 (14): 2382–83.

</div>

</div>

## Acknowledgements

This package was developed as part of an assessment for 2023 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `IntegraTRN` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.
