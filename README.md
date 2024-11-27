The files here contain raw data for each experiment in this paper and the R code analysis pipeline that was used to generate the final snRNA-seq dataset. Below is a brief description of each file.

**File 1: Zswim6_zenodo_NH.xlsx.** This is a spreadsheet containing the raw data for each graph in the manuscript. The data is organized by figure (see tabs at the bottom).

**File 2: individual 10X sample preprocessing_final.R.** After loading sparse gene expression matrices using the Read10X function, this R script was used to convert data from each biological 
sample into a Seurat object using the Seurat package (v4.3.0). In addition, the code normalizes and scales the raw gene expression data before performing PCA followed by Jackstraw analysis to 
determine the number of statistically significant PCs. 

**File 3: doublet_detection_allsamples.R.** This R script uses the doubletFinder_v3 package to detect putative doublets (see https://rdrr.io/github/chris-mcginnis-ucsf/DoubletFinder/ for details on package). 
Briefly, this R package generates artificial doublets based on gene expression data within the experimental data set, incorporates these artificial doublets into the real dataset, then predicts which experimental
cells/nuclei are doublets based on their proximity to the artificial doublets (i.e., a high proportion of neighboring artificial doublets would indicate a high probability that a cell/nucleus is a doublet).
This script is based on the workflow described in McGinnis et al., 2019 (reference below), using the conservative estimate of an 8% doublet rate. The values for pK, a parameter corresponding to the neighborhood
size in PC space used to calculate the number of artificial doublet nearest neighbors, are empirically determined for each data set. After processing all data sets and flagging putative doublets, the samples were merged
into a single seurat object which was run through the same doubletfinder pipeline (last 7 lines of code). This was done to ensure stringent criteria for removal of all doublets.



**DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors.**
McGinnis, Christopher S. et al.
Cell Systems, Volume 8, Issue 4, 329 - 337.e4
