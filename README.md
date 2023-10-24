
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IntegraTRN: Integrated omics for the inference of Transcriptional Regulatory Network

<!-- badges: start -->

[![License: GPL (\>=
3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://choosealicense.com/licenses/gpl-3.0/)
![GitHub issues](https://img.shields.io/github/issues/j-y26/IntegraTRN)
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
#> Loading required package: devtools
#> Loading required package: usethis
devtools::install_github("j-y26/IntegraTRN", build_vignettes = TRUE)
#> Downloading GitHub repo j-y26/IntegraTRN@HEAD
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\kirin\AppData\Local\Temp\Rtmp8Ctz4z\remotes6c686932563d\j-y26-IntegraTRN-e2ad6f2/DESCRIPTION' ...     checking for file 'C:\Users\kirin\AppData\Local\Temp\Rtmp8Ctz4z\remotes6c686932563d\j-y26-IntegraTRN-e2ad6f2/DESCRIPTION' ...   ✔  checking for file 'C:\Users\kirin\AppData\Local\Temp\Rtmp8Ctz4z\remotes6c686932563d\j-y26-IntegraTRN-e2ad6f2/DESCRIPTION' (567ms)
#>       ─  preparing 'IntegraTRN':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>       ─  building 'IntegraTRN_0.1.0.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/kirin/AppData/Local/Temp/RtmpGspdhQ/temp_libpath57106e6b626b'
#> (as 'lib' is unspecified)
library("IntegraTRN")
```

To run the shinyApp: `Under construction`

## Overview

``` r
ls("package:IntegraTRN")
#> [1] "protein_heart"             "protein_heart_samples"    
#> [3] "RNAseq_heart"              "RNAseq_heart_samples"     
#> [5] "smallRNAseq_heart"         "smallRNAseq_heart_samples"
data(package = "IntegraTRN")
browseVignettes("IntegraTRN")
#> No vignettes found by browseVignettes("IntegraTRN")
```

## Contributions

## References

## Acknowledgements

This package was developed as part of an assessment for 2023 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `IntegraTRN` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.
