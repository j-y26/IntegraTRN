# Purpose: Predictive estimation of regulatory small RNA - mRNA interaction
#          based on matched co-expression data
# Author: Jielin Yang
# Date: 2023-10-30
# Version: 1.0
# Bugs and Issues: None



















# Identify top genes that drive the differential expression of RNA or small RNAs
# top genes that drive the most variance separating the samples?

# use nearest neighbor matching to match the samples between the two groups
# then match the normalized counts of the top DE/variance genes between the two groups
# then use GENIE3 to predict the regulatory interactions between the top DE/variance genes and the small RNAs

# ensure when matching, the normalized counts of the small RNAs and the genes are generated from the same method