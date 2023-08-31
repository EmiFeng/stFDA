## Feature Differential Analysis for Spatial Transcriptomics (stFDA)

### 1. Introduction

stFDA is a computational tookit designed for the analysis of multiple spatial transcriptomics samples from 10X Visium. It provides methods to merge samples, extract cell type specific signatures from scRNA-seq data, spatial feature scoring, feature enriched spots detection and spatial differential analysis. 

### 2. Installation & Setup

To run through the analysis and plot results, simply download the R functions and load packages into your R session:

```R
Session Info:
R version 4.1.1
ggpubr_0.4.0
ggplot2_3.3.5
Seurat_4.0.5
SeuratObject_4.0.4 
irGSEA_1.1.2
```

```R
library(irGSEA)
library(Seurat)
library(ggplot2)
library(ggpubr)
source('stFDA_functions.R')
```

### 3. Run stFDA
```
# install packages from Bioconductor
bioconductor.packages <- c("GSEABase", "AUCell", "SummarizedExperiment", 
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree", 
                           "Nebulosa")
if (!requireNamespace(bioconductor.packages, quietly = TRUE)) { 
    BiocManager::install(bioconductor.packages, ask = F, update = F)
}

if (!requireNamespace("UCell", quietly = TRUE)) { 
    devtools::install_github("carmonalab/UCell")
}
if (!requireNamespace("irGSEA", quietly = TRUE)) { 
    devtools::install_github("chuiqin/irGSEA")
}
```
#### Merge datasets

```
rds <- c("S1.rds","S2.rds","S3.rds","S4.rds","S5.rds","S6.rds")
sample <- c('S1','S2','S3','S4','S5','S6')
group <- c('G1','G1','G1','G2','G2','G2')
PRO <- readRDS(rds[1])
obj_list <- list(PRO)
for (i in 2:length(rds)){
obj_list[i] <- readRDS(rds[i])
}
strds <- stFDA_mergeST(obj_list, sample_list = sample, group_list = group)
```

