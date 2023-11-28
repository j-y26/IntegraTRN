``` r
# Load the package
library(IntegraTRN)
library(dplyr)
```

## Introduction

Understanding the molecular mechanisms underlying complex diseases or
during the course of development has been a major goal in biomedical
research. The regulation of the transcriptome, both at the
transcriptional and post- transcriptional levels, is a key component of
these mechanisms. However, the complexity of the regulatory networks, in
particular as it involves multiple epigenetic players, renders the
dissection of these mechanisms challenging. In particular, the most
common may to construct a transcriptional regulatory network (TRN)
involves using known interactions at a certain condition. However, this
method is insufficient to explain the dynamics of the TRN in terms of
the changes observed between biological conditions. In this vignette, we
demonstrate how to use the **`IntegraTRN`** package to construct a TRN
using multi-omics data.

The package **`IntegraTRN`** is designed to integrate differential
analysis based on multiple omics data and construct a TRN based on both
curated and predicted interactions that underlie the differential
expression of genes. In the mean time, the curation of multi-omics data
could be difficult, so the package provides a flexibility to allow users
to decide what data can be integrated based on their own knowledge. The
analysis pipeline can be primarily divided into two parts: 1)
differential analysis of omics data; and 2) construction of the TRN. In
this vignette, we will demonstrate how to use the package to construct a
TRN using **transcriptomic**, **small RNAomic**, and **chromatin
accessibility** data.

See `help(package = "IntegraTRN")` for more information about the
package and `citation("IntegraTRN")` for how to cite the package. To
download the package, please use the following command:

``` r
require("devtools")
devtools::install_github("j-y26/IntegraTRN", build_vignettes = TRUE)
library("IntegraTRN")
```

To list all available functions in the package, use the following
command:

``` r
ls("package:IntegraTRN")
```

## Data preparation

The package has provided datasets for demonstration purposes. The
datasets are fully simulated or partially simulated based on real data
(Aharon-Yariv et al. 2023). Here is an overview of the basic omics data
used in this vignette:

``` r
data("RNAseq_heart")  # RNAseq count matrix
data("RNAseq_heart_samples")  # RNAseq sample information
data("smallRNAseq_heart")  # small RNAseq count matrix
data("smallRNAseq_heart_samples")  # small RNAseq sample information
```

These are RNAseq and small RNAseq data from human fetal heart tissues,
as well as the information about the samples.

Please use `?<dataset>` to see the details of each dataset.

We will now illustrate how to use the package to construct a TRN using
these datasets.

## Differential analysis of omics data

The first step of the analysis pipeline is to perform differential
analysis separately on each omic data. During the differential analysis,
we can retrieve insights into the “changes” that happen during the
biological process separately for each omic data. This exploratory
analysis is essential to understand what alterations in the
transcriptome can explain difference in the phenotype of interest, and
what variations in any of the regulatory players could potentially
contribute to these changes.

### Generate **`MOList`** object using multi omics data

To get started, we first need to generate a **`MOList`** object. The
**`MOList`** object stands for **`multi-omics list`**, which is the core
of entire analysis pipeline as it preserves our input data and the
results of the differential analysis in a centralized object. The
**`MOList`** object can be generated using the **`MOList`** function.

The **`MOList`** constructor has a large flexibility in terms of what
input omics data are available. At the minimum, when constructing a
**`MOList`** object for the first time, we need to provide the RNAseq
data and the sample grouping information. For any count-based omics data
(RNAseq, small RNAseq, and proteomics), the constructor takes a numeric
matrix of counts, where each row represents a gene or a protein, and
each column represents a sample. The grouping information can be any
vector of the same length as the number of samples, in the same order to
the columns of the count matrix. A continuous or a binary variable can
be used. However, we need to make sure we are supplying the same
grouping criteria for all omics data. Here we create a **`myMOList`**
using only RNAseq and small RNAseq data, using the **`age`** of the
samples as the grouping criteria for downstream differential analysis.

``` r
myMOList <- MOList(
  RNAseq = RNAseq_heart,
  RNAGroupBy = RNAseq_heart_samples$Age,
  smallRNAseq = smallRNAseq_heart,
  smallRNAGroupBy = smallRNAseq_heart_samples$Age
)
```

As long as a **`MOList`** object is constructed with RNAseq data,
additional omics data can be added to the existing object using the same
syntax. Any existing data can also be replaced by new data using the
same syntax.

Specifically for ATACseq data, we need to provide the path to the peak
files instead of the count matrix. Since most peak callers output peak
files in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html), we
can simply provide a path to a BED file, for each condition. Our package
has already provided two simulated peak files in the **`extdata`**
folder.

``` r
atacPeak1Path <- system.file("extdata", "peak1.bed", package = "IntegraTRN")
atacPeak2Path <- system.file("extdata", "peak2.bed", package = "IntegraTRN")
```

A consensus peak file should be provided for each condition. When
biological replicates are used for the ATACseq experiment, the peak
files should be merged into one and filtered for [blacklisted
regions](https://github.com/Boyle-Lab/Blacklist) before being provided
to the constructor (Amemiya, Kundaje, and Boyle 2019). Here, we simply
use the previously created **`myMOList`** object to append the ATACseq
data.

``` r
myMOList <- MOList(myMOList,
  pathATACpeak1 = atacPeak1Path,
  pathATACpeak2 = atacPeak2Path
)
```

### Differential analysis of omics data

Once the **`MOList`** object is created, we can perform differential
analysis on each omics data separately. However, since differential
analysis on count-based data is different that that of the epigenomic
data, the process is slightly different.

#### Differential analysis of count-based omics data

For count-based omics data, we can use the **`diffOmics`** function to
perform differential analysis. To simplify the process, the
**`diffOmics`** function essentially performs differential expression
analysis on all available count based omics data in the **`MOList`**
object. Since we can specify to use either **`DESeq2`** or **`edgeR`**
for the differential analysis, our underlying assumption is that the
expression data follows a negative binomial distribution (Love, Huber,
and Anders 2014) (Robinson, McCarthy, and Smyth 2010). Furthermore, the
differential expression analysis utilizes a generalized linear model
(GLM) to allow adjustment on covariates. In particular, the function
allows adjustment on the batch effect, which is a common issue in
next-generation sequencing (NGS) data.

``` r
myMOList <- diffOmics(myMOList,
  rnaseqBatch = RNAseq_heart_samples$Batch,
  program = "DESeq2"
)
#> Performing differential analysis on RNAseq data
#> 
#> Performing differential analysis on smallRNAseq data
```

The **`batch`** argument follow the same criteria as the grouping
information provided during the construction of the **`MOList`** object.
Therefore, this argument is not limited to the batch effect, but can be
any covariate that we want to adjust for.

``` r
# NOT run, just illustrate the usage of the batch argument
myMOList <- diffOmics(myMOList,
  smallRnaBatch = smallRNAseq_heart_samples$Sex,
  program = "edgeR"
)
```

##### Exploring differential gene expression analysis

To explore the results of the differential analysis, several
visualization methods are provided. First, to see the distribution of
differentially expressed mRNAs as well as exploring a optimal cutoff for
defining genes as differentially expressed, we can visualize the RNAseq
data using a volcano plot. We can specify the cutoffs for the adjusted
p-value and the log2 fold change (log2FC) to define differentially
expressed genes. Here, we define `log2FC > 0.1` and
`adjusted p-value < 0.05`, which is the default, as the cutoffs.

``` r
plotVolcanoRNA(myMOList, log2FC = 0.1)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-10-1.png" alt="Figure 1. Volcano plot of RNAseq data"  />
<p class="caption">
Figure 1. Volcano plot of RNAseq data
</p>

The plot also directly indicates the number of up- and down-regulated
genes.

For any count-based differential expression analysis, the result is
stored in a **`DETag`** object as an list element of the **`MOList`**
object. The **`myMOList`** object can simply be subsetted as a regular
list.

``` r
myRNADETag <- myMOList$DERNAseq
```

A simple snap shot of the `DETag` object can be directly invoked.

``` r
myRNADETag
#> A DETag S4 object
#> Method for differential analysis:  DESeq2 
#> 
#> Number of genes:  5000 
#> Number of genes with adjusted p-value < 0.05:  3223 
#> 
#> A snap shot for differential analysis results:
#>            baseMean log2FoldChange       lfcSE        stat       pvalue
#> B4GALT3   958.13887  -0.0078640979 0.009959527 -0.78960552 4.297582e-01
#> SLC26A11  458.38254  -0.0019483551 0.009851067 -0.19778113 8.432163e-01
#> BEAN1     114.23147  -0.0006515378 0.020920562 -0.03114342 9.751552e-01
#> STAT3    2121.37501  -0.0343187743 0.008212107 -4.17904601 2.927345e-05
#> STEAP3     56.69544   0.0311278206 0.016709689  1.86286052 6.248189e-02
#> PIP5K1C  1552.90544  -0.0978374269 0.012828549 -7.62653887 2.411402e-14
#>                  padj
#> B4GALT3  4.955699e-01
#> SLC26A11 8.747057e-01
#> BEAN1    9.810414e-01
#> STAT3    8.402252e-05
#> STEAP3   9.057972e-02
#> PIP5K1C  3.367881e-13
#> 
#> To access the full results, use the exportDE S4 method
```

Or if we want to retrieve the entire differential expression result, we
can use the **`exportDE`** function.

``` r
exportDE(myRNADETag) %>%
    head(5)  # For the first 5 genes
#>                  logFC       pvalue         padj
#> B4GALT3  -0.0078640979 4.297582e-01 4.955699e-01
#> SLC26A11 -0.0019483551 8.432163e-01 8.747057e-01
#> BEAN1    -0.0006515378 9.751552e-01 9.810414e-01
#> STAT3    -0.0343187743 2.927345e-05 8.402252e-05
#> STEAP3    0.0311278206 6.248189e-02 9.057972e-02
```

For publication purposes, we can also label some genes that we are
interested in on the volcano plot. Let’s first retrieve a few
differentially expressed genes

``` r
genesToLabel <- exportDE(myRNADETag) %>%
  dplyr::filter(padj < 0.05, abs(logFC) > 0.1) %>%
  head(3) %>%
  rownames()
```

Then we again plot the volcano plot, with some genes labeled.

We can make use of the **`ggrepel`** package to label the genes on the
plot.

``` r
rnaVPlot <- plotVolcanoRNA(myMOList, log2FC = 0.1)
rnaVPlot <- rnaVPlot + ggrepel::geom_label_repel(
  data = rnaVPlot$data[rownames(rnaVPlot$data) %in% genesToLabel, ],
  mapping = ggplot2::aes(
    x = logFC,
    y = -log10(padj),
    label = rownames(rnaVPlot$data[rownames(rnaVPlot$data) %in% genesToLabel, ])
  ),
  box.padding = 0.5,
  size = 3
)
rnaVPlot
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-15-1.png" alt="Figure 2. Volcano plot of RNAseq data with some genes labelled"  />
<p class="caption">
Figure 2. Volcano plot of RNAseq data with some genes labelled
</p>

##### Exploring differential small RNA expression analysis

Before moving on to small RNAseq data, we first need to annotate the
small RNAs for their types. To do so, we need to provide a data frame of
small RNA annotation, as in the name of the transcripts for each type of
small RNA. The following is an example:

``` r
smallRNAAnnotation <- data.frame(
  transcript = c("hsa-miR-1-3p", "hsa-miR-1-5p", "hsa-miR-2-5p"),
  type = c("miRNA", "miRNA", "miRNA")
)
```

To annotate the small RNAseq data, we can use the **`annotateSmallRNA`**
function.

``` r
# NOT run, just illustrate the usage of the smallRNAAnnotation argument
myMOList <- annotateSmallRNA(myMOList, anno = smallRNAAnnotation)
```

We need to make sure that the names of the transcripts are consistent
with the names in the count matrix, since different annotation databases
may provide different names for the same transcript.

Alternatively, the package internally provides an annotation database
for human small RNAs. This database is curated from several small RNA
annotation databases, including miRBase, circBase, piRBase, piRNABank,
GtRNAdb, GENCODE, and piRNACluster. Since we are dealing with human
data, we can make use of this internal database.

``` r
myMOList <- annotateSmallRNA(myMOList, anno = "human")
```

The **`annotateSmallRNA`** function will automatically checks for the
coverage of each small RNA type and will alert the user if the
annotation does not cover all existing transcripts.

Similar to the RNAseq data, we can visualize the distribution of
differentially expressed small RNAs using a volcano plot. At this point,
each type of small RNAs is colored separately.

``` r
# Plot Volcano plot of small RNAseq data
plotVolcanoSmallRNA(myMOList, log2FC = 0)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-19-1.png" alt="Figure 3. Volcano plot of small RNAseq data"  />
<p class="caption">
Figure 3. Volcano plot of small RNAseq data
</p>

Furthermore, although the function of microRNAs (miRNAs) is well
studied, the understanding of other types of small RNAs such as
piwi-interacting RNAs (piRNAs) and small nucleolar RNAs (snoRNAs) is
still limited. Therefore, we can assess the contribution of each type of
small RNA to the differences observed between samples. To do so, we can
perform a principal component analysis (PCA) on each type of small RNA
separately.

``` r
pcaPlotList <- plotSmallRNAPCAs(myMOList)
```

The **`plotSmallRNAPCAs`** function returns a list of PCA plots, one for
each type of small RNA. We can check the length of the list to see how
many PCA plots are generated.

``` r
length(pcaPlotList)
#> [1] 6
```

In this case, we have 6 PCA plots, one for each small RNA type. We can
plot the PCA plot for miRNA using the following command:

``` r
pcaPlotList$miRNA
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-22-1.png" alt="Figure 4. PCA plot of miRNA"  />
<p class="caption">
Figure 4. PCA plot of miRNA
</p>

Or, we can generate a combined PCA plot for all small RNA types by
making use of the **`ggpubr`** package.

``` r
ggpubr::ggarrange(pcaPlotList$miRNA, pcaPlotList$circRNA, pcaPlotList$piRNA, pcaPlotList$snoRNA,
    pcaPlotList$snRNA, pcaPlotList$tRNA, ncol = 3, nrow = 2)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-23-1.png" alt="Figure 5. Combined PCA plot of all small RNA types"  />
<p class="caption">
Figure 5. Combined PCA plot of all small RNA types
</p>

#### Differential analysis of chromatin accessibility

Analysis of epigenomic data requires several more steps than that of the
count-based omics data. In the case of ATACseq data, genomic regions are
profiled for their accessibility to the hyperactive Tn5 transposase.
This accessibility corresponds to the chromatin state of the region,
which can be either open or closed for DNA binding proteins,
particularly transcription factors (TFs). Therefore, to analyze
differential chromatin accessibility, we need to annotate the peaks with
the genes that are nearby. Furthermore, since we are interested in
transcriptional regulation, we can parse through the open chromatin
regions and identify the TF binding motifs that are enriched in these
regions. Differentially enriched motifs can then be used to infer the
TFs that are potentially involved in the regulation of the
differentially expressed genes.

To do so, we first need additional annotation information. Here, we make
use of the annotation database from the
**`TxDb.Hsapiens.UCSC.hg38.knownGene`** package, which is a **`TxDb`**
object that contains the genomic coordinates of all known genes in the
`hg38` genome assembly (Team and Maintainer 2023). We also need to
provide a more general annotation database, which is the
**`org.Hs.eg.db`** package to extract the gene information (Carlson
2023). Finally, scanning the open chromatin regions for TF binding
motifs requires identifying the sequence of the open chromatin regions.
Here, we use the corresponding genome FASTA database from the
**`BSgenome.Hsapiens.UCSC.hg38`** package (Pagès 2023).

``` r
# Annotate ATAC Peaks
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("BSgenome.Hsapiens.UCSC.hg38")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
```

Motif enrichment analysis would also require a known set of position
weight matrices (PWMs) for all known curated TFs. The PWM is a matrix
that contains the probability of each nucleotide (A, T, C, G) at every
position of the binding motif. The **`IntegraTRN`** package has already
provided a list of PWMs from the **`JASPAR2022`** database for
vertebrates. We can load the PWMs using the following command:

``` r
# Load PWMs from JASPAR
data("jasparVertebratePWM")
```

This is a **`PWMatrixList`** object defined in the **`TFBSTools`**
package (Tan and Lenhard 2016). This package also provides additional
functionality to retrieve PWMs from other databases or tools, such as
HOMER etc. Here, we will use our package-provided PWMs from
**`JASPAR2022`** (Castro-Mondragon et al. 2022).

The **`annotateATACPeaksMotif`** function performs both peak annotation
and motif enrichment analysis via a single call. Here, we define the
promoter region as the region **3kb** upstream and **3kb** downstream of
the transcription start site.

``` r
myMOList <- annotateATACPeaksMotif(myMOList,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = annoDb,
  bsgenome = bsgenome,
  pwmL = jasparVertebratePWM
)
```

##### Exploring differential chromatin accessibility

After performing peak annotation and motif enrichment analysis, we can
explore the analysis via several methods.

First, let us take a look at the annotation on the peaks for their
genomic features. This pie chart shows the percentage of peaks that are
annotated as each genomic feature.

``` r
plotATACAnno(myMOList)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-27-1.png" alt="Figure 6. Annotation of ATAC peaks"  />
<p class="caption">
Figure 6. Annotation of ATAC peaks
</p>

To further access the quality of the ATACseq experiment, we can plot the
coverage of the ATACseq peaks on the genome. This plot shows the where
the peaks are located in each chromosome.

``` r
plotATACCoverage(myMOList)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-28-1.png" alt="Figure 7. ATACseq coverage"  />
<p class="caption">
Figure 7. ATACseq coverage
</p>

On the left of the plots shows all differentially enriched motifs.

Then we can explore the key part of the analysis, which is the
differentially enriched motifs. Typically, the differential motif
enrichment can be very stringent especially when the ATACseq libraries
are not sequenced very deeply (less than 100 million reads) (Yan et al.
2020). Therefore, we can use a more relaxed cutoff for the p-value. The
**`plotATACMotifHeatmap`** function by default utilizes adjusted
p-value, but when a unadjusted p-value is provided, it will use the
unadjusted one and ignore the default. Since our ATACseq data is
simulated, here we use a relaxed cutoff.

``` r
plotATACMotifHeatmap(myMOList, pValue = 0.01)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-29-1.png" alt="Figure 8. Motif enrichment on ATACseq peaks"  />
<p class="caption">
Figure 8. Motif enrichment on ATACseq peaks
</p>

## Constructing the TRN

Exploring the differential analysis results can provide insights into
the changes that happen during the biological process. In addition, by
exploring these results, have a better understanding of the data and the
quality of the experiment. These understanding should allow us to make
better decisions on deciding the cutoffs used for selecting genes and
transcripts for constructing the TRN.

Now we can move on to the second part of the analysis pipeline, which is
to construct the TRN.

### Match samples between omics data

However, before directly construct the TRN, some data preprocessing is
required. First, the package analyzes predicted interactions by using a
tree-based algorithm to effectively predict regulatory interactions
based on co-expression between genes and small RNAs (Huynh-Thu et al.
2010). However, this algorithm requires that the samples between the
RNAseq and small RNAseq data are matched so that co-expression profiles
can be extracted. Therefore, we need to match the samples between the
RNAseq and small RNAseq data (Stuart et al. 2011).

We can perform an one-to-one optimal matching between the samples using
the `matchSamplesRNAsmallRNA` function. When we have further information
about the samples for RNAseq and small RNAseq, we can supply the
information as data frames. The algorithm by default extracts all common
variables between the two data frames and use them for matching.
Alternatively, we can also specify the variables used for matching using
the `varMatch` argument.

Here, we use all common variables, which are `Age` and `Sex` for
matching.

``` r
myMOList <- matchSamplesRNAsmallRNA(myMOList,
  sampleDFRNAseq = RNAseq_heart_samples,
  sampleDFSmallRNAseq = smallRNAseq_heart_samples
)
```

Let us take a look at the matched samples.

``` r
exportMatchResult(myMOList) %>%
    head(5)  # For the first 5 matches
#>      RNAseq SmallRNAseq
#> 1 sample_48    sample_1
#> 2 sample_51    sample_2
#> 3  sample_5    sample_3
#> 4 sample_49    sample_4
#> 5 sample_10    sample_5
```

There are a total of `37` matches.

``` r
nrow(exportMatchResult(myMOList))
#> [1] 37
```

### Using externally curated interactions

Depending on the biological question, we may want to use externally
curated interactions to be part of the TRN. For example, the miRNet
tools provides a comprehensive database of miRNA-target interactions,
but the information are usually not tissue-specific or all relevant to
the biological question (Chang et al. 2020). However, **`IntegraTRN`**
provides a flexibility to allow users to import externally curated
interactions, and during the construction of the TRN, these external
interactions are intersected with the predicted interactions as well as
the differential genes/enriched motifs to provide a more focused TRN.

Let us load some simulated miRNA-target interactions and TF-target
interactions. These data has been provided as part of the package. These
interactions are in the adjacency list format, where it is a list of two
vectors, regulator and target.

``` r
data("miR2Genes")
data("tf2Genes")
myMOList <- loadExtInteractions(myMOList, miR2Genes = miR2Genes, tf2Genes = tf2Genes)
```

### Setting cutoffs for omics data

As mentioned previously, the cutoffs used for selecting genes and
transcripts for constructing the TRN are important. Due to the large
number of omics data, a separate function is provide for users to set
the cutoffs. Here we use the cutoffs defined earlier in the differential
analysis.

``` r
omiCutoffs <- setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0.1,
  rnaTopGenes = 50,
  smallRNAAdjPval = 0.05,
  smallRNALogFC = 0.1,
  smallRNATopGenes = 50,
  atacMotifPval = 0.01
)
```

This can also be easily checked.

``` r
omiCutoffs
#> RNAseq adjusted p-value cutoff: < 0.05 
#> RNAseq log fold change cutoff: > 0.1 
#> RNAseq selecting top 50 DE genes 
#> smallRNAseq adjusted p-value cutoff: < 0.05 
#> smallRNAseq log fold change cutoff: > 0.1 
#> smallRNAseq selecting top 50 DE genes 
#> Proteomics adjusted p-value cutoff: < 0.05 
#> Proteomics log fold change cutoff: > 1 
#> ATACseq motif p-value cutoff: < 0.01 
#> ATACseq motif log fold enrichment cutoff: > 0
```

### Constructing and exploring the TRN

Finally, we can construct the TRN using the **`constructTRN`** function.
This function contains the logic to deal with different availability of
the omics data as well as the availability of external interactions. We
can specify whether we want to focus only on the up or down regulated
genes, or both. We also need to specify whether we want to use the
predicted interactions or not.

Here, we construct a network using both predicted and curated
interactions, for all differentially expressed genes that are both up-
and down-regulated. When predicted interactions are used, we can make
the process faster by using multiple threads, and a random seed can be
specified or reproducibility (Huynh-Thu et al. 2010).

``` r
myTRNet <- constructTRN(myMOList,
  omiCutoffs = omiCutoffs,
  smallRNAtypes = "all",
  targetDirection = "both",
  predicted = TRUE,
  nthreads = 1,
  seed = 123
)
#> Constructing a combined expression matrix for GENIE3...
#> Criteria for defining the top differentially expressed genes:
#>   log2 fold change cutoff:  0.1 
#>   adjusted p-value cutoff:  0.05 
#>   top 50  genes selected
#> Criteria for defining the top differentially expressed small RNAs:
#>   log2 fold change cutoff:  0.1 
#>   adjusted p-value cutoff:  0.05 
#>   top 50  small RNAs selected
#> Running GENIE3 for predicted interactions...
#>   Number of trees:  1000 
#>   Number of threads:  1 
#>   Tree method:  RF 
#> This may take a while...
```

We can easily plot the network using the **`plotNetwork`** function. We
can decide whether we want an interactive network. For large networks,
the interactive network provides a better way to explore the
interactions (Allaire et al. 2017).

``` r
plotNetwork(myTRNet, interactive = TRUE)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-37-1.png" alt="Figure 9. The Transcriptional Regulatory Network"  />
<p class="caption">
Figure 9. The Transcriptional Regulatory Network
</p>

When hovering over the nodes, we can see the name of the element. The
each type of the element (target gene, miRNA, piRNA, TF etc.) are color
coded. The layout of the network is computed by a force-directed
algorithm, but we can simply drag the nodes to better visualize the
interactions.

We can also plot a static network using the **`plotNetwork`** function.
But this time, we only plot the up-regulated genes and miRNAs with a
more stringent cutoff.

``` r
omiCutoffs2 <- omiCutoffs
omiCutoffs2$rnaAdjPval <- 0.01
omiCutoffs2$rnaLogFC <- 0.2
```

Let us construct a new TRN.

``` r
myTRNet2 <- constructTRN(myMOList,
  omiCutoffs = omiCutoffs2,
  smallRNAtypes = "miRNA",
  targetDirection = "up",
  predicted = TRUE
)
#> Constructing a combined expression matrix for GENIE3...
#> Criteria for defining the top differentially expressed genes:
#>   log2 fold change cutoff:  0.2 
#>   adjusted p-value cutoff:  0.01 
#>   top 50  genes selected
#> Criteria for defining the top differentially expressed small RNAs:
#>   log2 fold change cutoff:  0.1 
#>   adjusted p-value cutoff:  0.05 
#>   top 50  small RNAs selected
#> Running GENIE3 for predicted interactions...
#>   Number of trees:  1000 
#>   Number of threads:  1 
#>   Tree method:  RF 
#> This may take a while...
```

Then we can plot the network as a static figure (Csardi, Nepusz, et al.
2006).

``` r
plotNetwork(myTRNet2, vertex.size = 10)
```

<img src="F:/IntegraTRN/vignettes/Integrating_multi-omics_for_constructing_transcription_regulatory_networks_files/figure-markdown_github/unnamed-chunk-40-1.png" alt="Figure 10. The Transcriptional Regulatory Network of Up-regulated Genes"  />
<p class="caption">
Figure 10. The Transcriptional Regulatory Network of Up-regulated Genes
</p>

Some times it is hard to use the default parameters to plot a
publication-level network. However, users can customize the network
using the **`igraph`** package. We can extract the network as an
**`igraph`** object and performs further customization (Csardi, Nepusz,
et al. 2006).

``` r
exportIgraph(myTRNet2)
#> IGRAPH cc3b273 DN-B 73 141 -- 
#> + attr: name (v/c), type (v/c)
#> + edges from cc3b273 (vertex names):
#>  [1] hsa-miR-1301-3p ->KLHL40  hsa-mir-183     ->KBTBD12
#>  [3] hsa-miR-323a-3p ->TSPO    hsa-miR-1301-3p ->LRRC14B
#>  [5] hsa-miR-6747-3p ->COX7C   hsa-miR-1910-5p ->NRAP   
#>  [7] hsa-miR-31-5p   ->COX7C   hsa-miR-323a-3p ->DBP    
#>  [9] hsa-miR-323a-3p ->C1QL1   hsa-miR-1301-3p ->CEBPD  
#> [11] hsa-miR-31-5p   ->CEBPD   hsa-miR-323a-3p ->FAU    
#> [13] hsa-miR-323a-3p ->G0S2    hsa-mir-130b    ->LGALS3 
#> [15] hsa-mir-183     ->FAM129A hsa-miR-6747-3p ->CEBPD  
#> + ... omitted several edges
```

## References

Aharon-Yariv, Adar, Yaxu Wang, Abdalla Ahmed, and Paul Delgado-Olguı́n.
2023. “Integrated Small RNA, mRNA and Protein Omics Reveal a miRNA
Network Orchestrating Metabolic Maturation of the Developing Human
Heart.” *BMC Genomics* 24 (1): 1–18.

Allaire, JJ, Peter Ellis, Christopher Gandrud, Kevin Kuo, BW Lewis,
Jonathan Owen, Kenton Russell, et al. 2017. *networkD3: D3 JavaScript
Network Graphs from r*. <https://CRAN.R-project.org/package=networkD3>.

Amemiya, Haley M, Anshul Kundaje, and Alan P Boyle. 2019. “The ENCODE
Blacklist: Identification of Problematic Regions of the Genome.”
*Scientific Reports* 9 (1): 9354.

Carlson, Marc. 2023. *Org.hs.eg.db: Genome Wide Annotation for Human*.

Castro-Mondragon, Jaime A, Rafael Riudavets-Puig, Ieva Rauluseviciute,
Roza Berhanu Lemma, Laura Turchi, Romain Blanc-Mathieu, Jeremy Lucas, et
al. 2022. “JASPAR 2022: The 9th Release of the Open-Access Database of
Transcription Factor Binding Profiles.” *Nucleic Acids Research* 50
(D1): D165–73.

Chang, Le, Guangyan Zhou, Othman Soufan, and Jianguo Xia. 2020. “miRNet
2.0: Network-Based Visual Analytics for miRNA Functional Analysis and
Systems Biology.” *Nucleic Acids Research* 48 (W1): W244–51.

Csardi, Gabor, Tamas Nepusz, et al. 2006. “The Igraph Software Package
for Complex Network Research.” *InterJournal, Complex Systems* 1695 (5):
1–9.

Huynh-Thu, Vân Anh, Alexandre Irrthum, Louis Wehenkel, and Pierre
Geurts. 2010. “Inferring Regulatory Networks from Expression Data Using
Tree-Based Methods.” *PloS One* 5 (9): e12776.

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 1–21.

Pagès, Hervé. 2023. *BSgenome: Software Infrastructure for Efficient
Representation of Full Genomes and Their SNPs*.
<https://doi.org/10.18129/B9.bioc.BSgenome>.

Robinson, Mark D, Davis J McCarthy, and Gordon K Smyth. 2010. “edgeR: A
Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data.” *Bioinformatics* 26 (1): 139–40.

Stuart, Elizabeth A, Gary King, Kosuke Imai, and Daniel Ho. 2011.
“MatchIt: Nonparametric Preprocessing for Parametric Causal Inference.”
*Journal of Statistical Software*.

Tan, Ge, and Boris Lenhard. 2016. “TFBSTools: An r/Bioconductor Package
for Transcription Factor Binding Site Analysis.” *Bioinformatics* 32
(10): 1555–56.

Team, Bioconductor Core, and Bioconductor Package Maintainer. 2023.
*TxDb.hsapiens.UCSC.hg38.knownGene: Annotation Package for TxDb
Object(s)*.

Villanueva, Randle Aaron M, and Zhuo Job Chen. 2019. “Ggplot2: Elegant
Graphics for Data Analysis.” Taylor & Francis.

Wickham, Hadley, Romain François, Lionel Henry, Kirill Müller, and Davis
Vaughan. 2023. *Dplyr: A Grammar of Data Manipulation*.
<https://CRAN.R-project.org/package=dplyr>.

Yan, Feng, David R Powell, David J Curtis, and Nicholas C Wong. 2020.
“From Reads to Insight: A Hitchhiker’s Guide to ATAC-Seq Data Analysis.”
*Genome Biology* 21: 1–16.

## Session information

``` r
sessionInfo("IntegraTRN")
#> R version 4.3.1 (2023-06-16 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 11 x64 (build 22621)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=C                    LC_CTYPE=English_Canada.utf8   
#> [3] LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C                   
#> [5] LC_TIME=English_Canada.utf8    
#> 
#> time zone: America/Toronto
#> tzcode source: internal
#> 
#> attached base packages:
#> character(0)
#> 
#> other attached packages:
#> [1] IntegraTRN_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] IRanges_2.36.0                          
#>   [2] R.methodsS3_1.8.2                       
#>   [3] progress_1.2.2                          
#>   [4] vsn_3.70.0                              
#>   [5] urlchecker_1.0.1                        
#>   [6] poweRlaw_0.70.6                         
#>   [7] Biostrings_2.70.1                       
#>   [8] TxDb.Hsapiens.UCSC.hg38.knownGene_3.18.0
#>   [9] vctrs_0.6.4                             
#>  [10] digest_0.6.33                           
#>  [11] png_0.1-8                               
#>  [12] shape_1.4.6                             
#>  [13] ggrepel_0.9.4                           
#>  [14] MASS_7.3-60                             
#>  [15] reshape2_1.4.4                          
#>  [16] httpuv_1.6.11                           
#>  [17] foreach_1.5.2                           
#>  [18] BiocGenerics_0.48.1                     
#>  [19] qvalue_2.34.0                           
#>  [20] withr_2.5.2                             
#>  [21] stabs_0.6-4                             
#>  [22] xfun_0.40                               
#>  [23] ggfun_0.1.3                             
#>  [24] ggpubr_0.6.0                            
#>  [25] ellipsis_0.3.2                          
#>  [26] survival_3.5-5                          
#>  [27] memoise_2.0.1                           
#>  [28] utils_4.3.1                             
#>  [29] profvis_0.3.8                           
#>  [30] tidytree_0.4.5                          
#>  [31] networkD3_0.4                           
#>  [32] zoo_1.8-12                              
#>  [33] GlobalOptions_0.1.2                     
#>  [34] gtools_3.9.4                            
#>  [35] DNAcopy_1.76.0                          
#>  [36] R.oo_1.25.0                             
#>  [37] prettyunits_1.2.0                       
#>  [38] KEGGREST_1.42.0                         
#>  [39] promises_1.2.1                          
#>  [40] httr_1.4.7                              
#>  [41] rstatix_0.7.2                           
#>  [42] restfulr_0.0.15                         
#>  [43] GENIE3_1.24.0                           
#>  [44] ps_1.7.5                                
#>  [45] rstudioapi_0.15.0                       
#>  [46] miniUI_0.1.1.1                          
#>  [47] generics_0.1.3                          
#>  [48] DOSE_3.28.0                             
#>  [49] processx_3.8.2                          
#>  [50] curl_5.1.0                              
#>  [51] S4Vectors_0.40.1                        
#>  [52] zlibbioc_1.48.0                         
#>  [53] ggraph_2.1.0                            
#>  [54] polyclip_1.10-6                         
#>  [55] GenomeInfoDbData_1.2.11                 
#>  [56] SparseArray_1.2.0                       
#>  [57] interactiveDisplayBase_1.40.0           
#>  [58] xtable_1.8-4                            
#>  [59] stringr_1.5.0                           
#>  [60] desc_1.4.2                              
#>  [61] pracma_2.4.2                            
#>  [62] doParallel_1.0.17                       
#>  [63] evaluate_0.23                           
#>  [64] S4Arrays_1.2.0                          
#>  [65] BiocFileCache_2.10.1                    
#>  [66] preprocessCore_1.64.0                   
#>  [67] hms_1.1.3                               
#>  [68] glmnet_4.1-8                            
#>  [69] GenomicRanges_1.54.1                    
#>  [70] colorspace_2.1-0                        
#>  [71] filelock_1.0.2                          
#>  [72] webshot2_0.1.1                          
#>  [73] graphics_4.3.1                          
#>  [74] magrittr_2.0.3                          
#>  [75] readr_2.1.4                             
#>  [76] later_1.3.1                             
#>  [77] viridis_0.6.4                           
#>  [78] ggtree_3.10.0                           
#>  [79] lattice_0.21-8                          
#>  [80] genefilter_1.84.0                       
#>  [81] methods_4.3.1                           
#>  [82] XML_3.99-0.15                           
#>  [83] shadowtext_0.1.2                        
#>  [84] cowplot_1.1.1                           
#>  [85] matrixStats_1.0.0                       
#>  [86] pillar_1.9.0                            
#>  [87] nlme_3.1-162                            
#>  [88] iterators_1.0.14                        
#>  [89] caTools_1.18.2                          
#>  [90] compiler_4.3.1                          
#>  [91] stringi_1.7.12                          
#>  [92] stats_4.3.1                             
#>  [93] SummarizedExperiment_1.32.0             
#>  [94] devtools_2.4.5                          
#>  [95] GenomicAlignments_1.38.0                
#>  [96] MPO.db_0.99.7                           
#>  [97] plyr_1.8.9                              
#>  [98] crayon_1.5.2                            
#>  [99] abind_1.4-5                             
#> [100] BiocIO_1.12.0                           
#> [101] truncnorm_1.0-9                         
#> [102] gridGraphics_0.5-1                      
#> [103] sm_2.2-5.7.1                            
#> [104] chk_0.9.1                               
#> [105] locfit_1.5-9.8                          
#> [106] graphlayouts_1.0.2                      
#> [107] org.Hs.eg.db_3.18.0                     
#> [108] bit_4.0.5                               
#> [109] chromote_0.1.2                          
#> [110] dplyr_1.1.3                             
#> [111] fastmatch_1.1-4                         
#> [112] codetools_0.2-19                        
#> [113] Ringo_1.66.0                            
#> [114] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 
#> [115] GetoptLong_1.0.5                        
#> [116] mime_0.12                               
#> [117] splines_4.3.1                           
#> [118] circlize_0.4.15                         
#> [119] Rcpp_1.0.11                             
#> [120] datasets_4.3.1                          
#> [121] dbplyr_2.4.0                            
#> [122] HDO.db_0.99.1                           
#> [123] knitr_1.45                              
#> [124] blob_1.2.4                              
#> [125] utf8_1.2.4                              
#> [126] clue_0.3-65                             
#> [127] BiocVersion_3.18.0                      
#> [128] seqLogo_1.68.0                          
#> [129] fs_1.6.3                                
#> [130] Rdpack_2.6.9000                         
#> [131] Repitools_1.48.0                        
#> [132] pkgbuild_1.4.2                          
#> [133] ggsignif_0.6.4                          
#> [134] ggplotify_0.1.2                         
#> [135] tibble_3.2.1                            
#> [136] Matrix_1.6-1.1                          
#> [137] callr_3.7.3                             
#> [138] statmod_1.5.0                           
#> [139] tzdb_0.4.0                              
#> [140] tweenr_2.0.2                            
#> [141] pkgconfig_2.0.3                         
#> [142] tools_4.3.1                             
#> [143] cachem_1.0.8                            
#> [144] rbibutils_2.2.16                        
#> [145] RSQLite_2.3.3                           
#> [146] viridisLite_0.4.2                       
#> [147] DBI_1.1.3                               
#> [148] monaLisa_1.8.0                          
#> [149] fastmap_1.1.1                           
#> [150] rmarkdown_2.25                          
#> [151] scales_1.2.1                            
#> [152] grid_4.3.1                              
#> [153] usethis_2.2.2                           
#> [154] Rsamtools_2.18.0                        
#> [155] broom_1.0.5                             
#> [156] AnnotationHub_3.10.0                    
#> [157] patchwork_1.1.3                         
#> [158] BiocManager_1.30.22                     
#> [159] carData_3.0-5                           
#> [160] farver_2.1.1                            
#> [161] tidygraph_1.2.3                         
#> [162] scatterpie_0.2.1                        
#> [163] yaml_2.3.7                              
#> [164] roxygen2_7.2.3                          
#> [165] MatrixGenerics_1.14.0                   
#> [166] rtracklayer_1.62.0                      
#> [167] cli_3.6.1                               
#> [168] purrr_1.0.2                             
#> [169] webshot_0.5.5                           
#> [170] stats4_4.3.1                            
#> [171] lifecycle_1.0.4                         
#> [172] Biobase_2.62.0                          
#> [173] sessioninfo_1.2.2                       
#> [174] backports_1.4.1                         
#> [175] BSgenome.Hsapiens.UCSC.hg38_1.4.5       
#> [176] BiocParallel_1.36.0                     
#> [177] annotate_1.80.0                         
#> [178] gtable_0.3.4                            
#> [179] rjson_0.2.21                            
#> [180] base_4.3.1                              
#> [181] ChIPseeker_1.38.0                       
#> [182] parallel_4.3.1                          
#> [183] ape_5.7-1                               
#> [184] limma_3.58.1                            
#> [185] jsonlite_1.8.7                          
#> [186] edgeR_4.0.1                             
#> [187] TFBSTools_1.40.0                        
#> [188] bitops_1.0-7                            
#> [189] ggplot2_3.4.4                           
#> [190] HPO.db_0.99.2                           
#> [191] bit64_4.0.5                             
#> [192] yulab.utils_0.1.0                       
#> [193] CNEr_1.38.0                             
#> [194] highr_0.10                              
#> [195] GOSemSim_2.28.0                         
#> [196] rlemon_0.2.1                            
#> [197] R.utils_2.12.2                          
#> [198] lazyeval_0.2.2                          
#> [199] shiny_1.7.5.1                           
#> [200] htmltools_0.5.6.1                       
#> [201] affy_1.80.0                             
#> [202] enrichplot_1.22.0                       
#> [203] GO.db_3.18.0                            
#> [204] rappdirs_0.3.3                          
#> [205] formatR_1.14                            
#> [206] glue_1.6.2                              
#> [207] TFMPvalue_0.0.9                         
#> [208] XVector_0.42.0                          
#> [209] MatchIt_4.5.5                           
#> [210] RCurl_1.98-1.13                         
#> [211] rprojroot_2.0.4                         
#> [212] treeio_1.26.0                           
#> [213] BSgenome_1.70.1                         
#> [214] Rsolnp_1.16                             
#> [215] gridExtra_2.3                           
#> [216] gsmoothr_0.1.7                          
#> [217] boot_1.3-28.1                           
#> [218] igraph_1.5.1                            
#> [219] R6_2.5.1                                
#> [220] tidyr_1.3.0                             
#> [221] DESeq2_1.42.0                           
#> [222] gplots_3.1.3                            
#> [223] labeling_0.4.3                          
#> [224] GenomicFeatures_1.54.1                  
#> [225] cluster_2.1.4                           
#> [226] pkgload_1.3.3                           
#> [227] aplot_0.2.2                             
#> [228] GenomeInfoDb_1.38.0                     
#> [229] vioplot_0.4.0                           
#> [230] DirichletMultinomial_1.44.0             
#> [231] DelayedArray_0.28.0                     
#> [232] tidyselect_1.2.0                        
#> [233] plotrix_3.8-3                           
#> [234] optmatch_0.10.6                         
#> [235] ggforce_0.4.1                           
#> [236] xml2_1.3.5                              
#> [237] car_3.1-2                               
#> [238] AnnotationDbi_1.64.1                    
#> [239] munsell_0.5.0                           
#> [240] KernSmooth_2.23-21                      
#> [241] affyio_1.72.0                           
#> [242] data.table_1.14.8                       
#> [243] grDevices_4.3.1                         
#> [244] websocket_1.4.1                         
#> [245] htmlwidgets_1.6.2                       
#> [246] fgsea_1.28.0                            
#> [247] ComplexHeatmap_2.18.0                   
#> [248] RColorBrewer_1.1-3                      
#> [249] biomaRt_2.58.0                          
#> [250] rlang_1.1.1                             
#> [251] remotes_2.4.2.1                         
#> [252] fansi_1.0.5
```
