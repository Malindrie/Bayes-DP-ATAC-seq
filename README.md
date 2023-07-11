# Bayes-DP-ATAC-seq
Code used for the single cell RNA-seq &amp; ATAC-seq analysis of the paper "Bayesian non parametric mixture models reveal modes of regulation in chromatin accessibility and identifies genes that define cell identity" by Dharmaratne et al. 

The Rscript [Pre-processing ATAC-seq](https://github.com/Malindrie/Bayes-DP-ATAC-seq/blob/main/Pre-processing%20ATAC-seq.R) outlines the pre-processing performed for RNA and ATAC data.

The Rscipt [Accessiblbe promoter regions](https://github.com/Malindrie/Bayes-DP-ATAC-seq/blob/main/Accessiblbe%20promoter%20regions.R) outlines how to identify the promoter peaks in the ATAC-seq data and assigning the nearest gene TSS.

In the Rscript [Creating pseudocells](https://github.com/Malindrie/Bayes-DP-ATAC-seq/blob/main/Creating%20pseudocells.R) outlines creation of pseudocells based on the RNA-seq profiles of the data.

The Rscript [Dirichlet process](https://github.com/Malindrie/Bayes-DP-ATAC-seq/blob/main/Dirichlet%20process.R) is where we have implemented Bayesian Dirichlet Processes for identifying the promoter modes. Here we also have implemented the diagnostic tests for identifyin the correct number of clusters.

The remaining Rcripts are for downstream analysis. 
