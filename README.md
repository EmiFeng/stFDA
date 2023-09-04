## Feature Differential Analysis for Spatial Transcriptomics (stFDA)

### 1. Introduction

stFDA is a computational tookit designed for the analysis of multiple spatial transcriptomics samples from 10X Visium. It provides methods to merge samples, extract cell type specific signatures from scRNA-seq data, spatial feature scoring, feature enriched spots detection and spatial differential analysis. 

### 2. Installation & Setup

To run through the analysis and plot results, simply download the R functions and load packages into your R session:

```R
Session Info:
R version 4.1.1
ggpubr
ggplot2
Seurat_4.0.5
SeuratObject_4.0.4 
irGSEA_1.1.2

# install packages from CRAN
cran.packages <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges", 
                   "gghalves", "Seurat", "SeuratObject", "methods", 
                   "devtools", "BiocManager","data.table","doParallel",
                   "doRNG","ggpubr","ggplot2")
if (!requireNamespace(cran.packages, quietly = TRUE)) { 
    install.packages(cran.packages, ask = F, update = F)
}

# install packages from Bioconductor
bioconductor.packages <- c("GSEABase", "AUCell", "SummarizedExperiment", 
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree", 
                           "Nebulosa","clusterProfiler")
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

```R
library(irGSEA)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)

```

### 3. Run stFDA

#### Quick start

functions :download the R function from github

Example data: download the Visium data of human white adipose tissues (hWATs) from published work (https://doi.org/10.1016/j.cmet.2021.07.018) and run Space Ranger & Seurat pipelines.

```
#load stFDA functions
source('stFDA_functions.R')
```

#### Merge datasets

Function: merge spatial Seurat objects

**obj_list:** list of st seurat objects. 

**sample_list:** sample names.

**group_list:** group names.

```R
sam_gr <- read.table('test/group_BMI.xls',header = T)

> sam_gr
   sample      group   strds
1     S41       lean S41.rds
2     S42 overweight S42.rds
3     S43 overweight S43.rds
4     S44       lean S44.rds
5     S45       lean S45.rds
6     S46 overweight S46.rds
7     S47 overweight S47.rds
8     S48       lean S48.rds
9     S49      obese S49.rds
10    S50      obese S50.rds
11    S51      obese S51.rds

strds <- stFDA_mergeST(obj_list = sam_gr$strds, sample_list = sam_gr$sample, group_list = sam_gr$group)
```

#### Cell type signature genesets

Function: to get cell type signature genesets from scRNA-seq data, based on Seurat FindMarkers() function. Return a list of signatures, can be used for function stFDA_Scoring()

**scrds:** scRNA-seq seurat object with cell type annotation. 

**number:** maximum number of genes selected for a signature, default = 50, which means select top 50 ( avg_log2FC ordered) upregulated DEGs with significance. Note that mitochondrial and ribosomal genes will be excluded for signature gene selection. 

**group.by:** whether to change active.ident, default not to change; if change, give the colname stored in meta.data, for example, 'new_ident'.

**downsample:** whether to downsample cells for idents, default not to; suggest to use when cell number is too large, for exsample, downsample = 2000.

```R
sc_signature <- stFDA_scSignatures(scrds, downsample = 500, number = 50) #downsample = 500 for test use.

save(sc_signature, file = 'sc_signature.Rdata') #save signature list as .Rdata file.

#load signature list file afterwards.
load('sc_signature.Rdata') 
```

```R
#signature of Endothelial Cells derived from scRNA-seq data
> sc_signature[['Endothelial cells']]
 [1] "ACTA2"    "PECAM1"   "AQP1"     "VWF"      "ACKR1"    "PLVAP"
 [7] "IGFBP7"   "TAGLN"    "A2M"      "NOTCH3"   "HLA-B"    "C2orf40"
[13] "HLA-DRB1" "ID1"      "RBP7"     "HLA-DRA"  "HLA-E"    "SPARCL1"
[19] "CLDN5"    "MT2A"     "IFI27"    "CD300LG"  "HLA-C"    "CD74"
[25] "EPAS1"    "RGS5"     "NAA38"    "ADGRF5"   "TINAGL1"  "BTNL9"
[31] "IFITM3"   "ID3"      "B2M"      "SOX18"    "BST2"     "VAMP5"
[37] "RTF2"     "CD93"     "ADGRL4"   "IFITM2"   "COL18A1"  "HLA-A"
[43] "ENG"      "FAU"      "PDGFRB"   "PSME1"    "EGFL7"    "TPM2"
[49] "RAMP2"    "HLA-DPA1"
```

#### Spatial Feature Scoring

stFDA_Scoring enables integrated spatial feature scoring. Feature scoring is based on irGSEA.score() function from irGSEA package, which collected 4 scoring methods, AUCell, Ucell, singscore and ssgsea. Returns a feature scoring table as well as spatial feature plots for features. Suggest not to input large amount of genesets for analysis.

**obj:** spatial seurat object, if you have multiple samples, you can use the result from stFDA_mergeST.

**features:** functional genesets, can be feature list (return from stFDA_scSignatures()) or gmt file. 

**method:** choose one scoring method from "AUCell", "UCell", "singscore", "ssgsea", default = "ssgsea".

**assay:** which assay to use, default = 'Spatial', check Seurat object for more information.

**slot:** which data slot to use, default = 'data', check Seurat object for more information.

**custom:** Default TRUE. Set it to TRUE when input own genesets; set it to FALSE when use msigdb genesets.

**msigdb:** Default FALSE. Set it to TRUE when custom = FALSE.

**category:** Choose a msigdb category for analysis, default='H' . Use msigdbr::msigdbr_collections to view all available collections gene sets.

**species:** Default Homo sapiens. Use msigdbr::msigdbr_show_species() to view all available species. The parameter works if msigdb is True.

**outdir:** output directory path.

**prefix:** prefix of output files.

```R
ssgsea_table <- stFDA_Scoring(obj = strds, features = sc_signature, gmt = FALSE, custom = TRUE, method = 'ssgsea',prefix = 'ssgsea', outdir = '.') #Custom analysis for cell type features extracted from stFDA_scSignatures().

ssgsea_table <- stFDA_Scoring(obj = strds, features = 'geneset.gmt', gmt = TRUE, custom = TRUE, method = 'ssgsea',prefix = 'ssgsea', outdir = '.') #Custom analysis for .gmt file. 

ssgsea_table <- stFDA_Scoring(obj = strds, custom = FALSE, method = 'ssgsea',msigdb = TRUE, category = 'H',prefix = 'ssgsea', species = 'Homo sapiens', outdir = '.') #Functional analysis of msigdb pathways.
```

```R
> ssgsea_table[1:5,1:3]
                       Mononuclear.phagocytes Adipocytes Endothelial.cells
S41_AAACACCAATAACTGC-1               3691.054   3782.827          3775.384
S41_AAACCGGGTAGGTACC-1               2144.613   2729.821          2261.770
S41_AAACCGTTCGTCCAGG-1               2529.504   3855.121          2915.817
S41_AAACCTCATGAAGTTG-1               2561.356   3876.238          2549.837
S41_AAACGAGACGGTTGAT-1               2342.973   3084.321          2202.788
```
![image](https://github.com/EmiFeng/stFDA/assets/41361743/358eaa1d-4def-4fba-9317-acee84a095f4)

#### Spatially Enriched Spots

Based on the table of feature scoring from stFDA_Scoring(), stFDA_threshold() and stFDA_enrichedSpots() could detect spatially enriched spots. Return assigned table as well as the fraction of enriched spots of every sample. For every celltype or feature (column), 1 represents for celltype/feature-enriched spots, while 0 represents for non-enriched spots. Multiple celltypes or features can be enriched in a single spot.

**data:** the table of feature scoring from stFDA_Scoring()

**obj:** spatial seurat object, if you have multiple samples, you can use the result from stFDA_mergeST.

**outlier:** to define a high-threshold as the maximum score after excluding top X percent outliers with highest scores. default = 0.01

**percent:** to determine a low-threshold cross which was regarded as a enriched spot. default = 0.7

**outdir:** output directory path.

**prefix:** prefix of output files.

```R
ssgsea_result <- stFDA_threshold(data = ssgsea_table, obj = strds, outlier = 0.01, percent = 0.7, outdir = ".")

> ssgsea_result[1:5,1:3]
                       Mononuclear phagocytes Adipocytes Endothelial cells
S41_AAACACCAATAACTGC-1                      1          0                 1
S41_AAACCGGGTAGGTACC-1                      0          0                 0
S41_AAACCGTTCGTCCAGG-1                      0          0                 0
S41_AAACCTCATGAAGTTG-1                      0          0                 0
S41_AAACGAGACGGTTGAT-1                      0          0                 0
```
![image](https://github.com/EmiFeng/stFDA/assets/41361743/34bc4223-a2f8-485c-9cde-e1968750eb9f)

```R
ssgsea_enriched <-  stFDA_enrichedSpots(assign = ssgsea_result, strds = strds, outdir = ".", prefix = 'ssgsea')

> head(ssgsea_enriched)
  sample      group   percent               celltype
1    S41       lean 0.1828428 Mononuclear phagocytes
2    S42 overweight 0.2011923 Mononuclear phagocytes
3    S43 overweight 0.2949355 Mononuclear phagocytes
4    S44       lean 0.2788967 Mononuclear phagocytes
5    S45       lean 0.3268136 Mononuclear phagocytes
6    S46 overweight 0.1671996 Mononuclear phagocytes
```

#### sample/group plot and statistical test

Draw sample barplot and group boxplot.  P-values from **t.test** or **wilcox.test** can be chosen for estimating the statistical significance.

**stat_test**: whether to do statistical test, choose TRUE when have biological repeats. default = FALSE

test: **t.test** or **wilcox.test**,  default = t.test

```R
stFDA_diffplot(table = ssgsea_enriched, outdir = '.', prefix = 'ssgsea', stat_test = TRUE, test = 't.test')

> sig
                     celltype     .y.     group1     group2           p p.adj
4                  Adipocytes percent       lean overweight 0.500492926 1.000
5                  Adipocytes percent       lean      obese 0.068152912 0.750
6                  Adipocytes percent overweight      obese 0.002379400 0.036
10 Adipose-derived stem cells percent       lean overweight 0.073270018 0.750
11 Adipose-derived stem cells percent       lean      obese 0.186254703 1.000
12 Adipose-derived stem cells percent overweight      obese 0.009148435 0.130
   p.format p.signif method
4    0.5005       ns T-test
5    0.0682       ns T-test
6    0.0024       ** T-test
10   0.0733       ns T-test
11   0.1863       ns T-test
12   0.0091       ** T-test
```
![image](https://github.com/EmiFeng/stFDA/assets/41361743/d19ced80-5cba-4b52-94e3-4056fb6c9b1c)


