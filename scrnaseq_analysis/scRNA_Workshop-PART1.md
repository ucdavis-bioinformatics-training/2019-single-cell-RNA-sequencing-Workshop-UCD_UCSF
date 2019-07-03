---
title: "Single Cell RNAseq Part 1"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

# Single Cell Analysis with Seurat and some custom code!

[Seurat](http://satijalab.org/seurat/) is a popular R package that is designed for QC, analysis, and exploration of single cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single cell transcriptomic measurements, and to integrate diverse types of single cell data. Further, the authors provide several [tutorials](http://satijalab.org/seurat/get_started.html), on their website.

Download and expand the expression_tables_cellrangerV3.zip file to extract the single cell matrix files for the three samples. These are isolated mouse cells ran on the 10X genomics platform (3' expression V2) for single cell RNA sequencing, sequenced with UC Davis on 1 HiSeq 4000 lane. The Experiment contains 3 samples, each merged from 2 original samples and then randomly subsamples to 1000 cells each.

The three samples are, Dorsal root ganglion neurons :
At weaning, Ttpa+/+ mice were fed a normal diet (35 mg of dl-α-tocopheryl acetate/kg, vitE+) while and Ttpa-/- mice were fed either an α-TOH-deficient diet (<10 mg of dl-α-tocopheryl acetate/kg, vitE-), or α-TOH-supplemented diet (600 mg of dl-α-tocopheryl acetate/kg, vitE+++) diet.

* UCD_Adj_VitE - normal + Vitamin E
* UCD_Supp_VitE - Vitamin E supplimented by diet.
* UCD_VitE_Def - Vitamin E deficient animals

We start with loading needed libraries for R, at this time all we need is the package [Seurat](http://satijalab.org/seurat/).

```r
library(Seurat)
```

## Load the Cell Ranger Matrix Data and create the base Seurat object.
Cell Ranger provides a function `cellranger aggr` that will combine multiple samples into a single matrix file. However, when processing data in R and Seurat this is unnecessary and we can aggregate them in R.

Seurat provides a function `Read10X` to read in 10X data folder. First we read in data from each individual sample folder. First, we initialize the Seurat object (`CreateSeuratObject`) with the raw (non-normalized data). Keep all genes expressed in >= 10 cells. Keep all cells with at least 200 detected genes. Also extracting sample names, calculating and adding in the metadata mitochondrial percentage of each cell. Adding in the metadata batchid and cell cycle. Finally, saving the raw Seurat object.


```r
dataset_loc <- "./expression_tables_cellrangerV3"
ids <- c("UCD_Adj_VitE", "UCD_Supp_VitE", "UCD_VitE_Def")

d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "scRNA workshop",
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "\\-")
```

### The percentage of reads that map to the mitochondrial genome

* Low-quality / dying cells often exhibit extensive mitochondrial contamination.
* We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features.
* We use the set of all genes, in mouse these genes can be identified as those that begin with 'mt', in human data they begin with MT.


```r
experiment.aggregate$percent.mito <- PercentageFeatureSet(experiment.aggregate, pattern = "^mt-")
```


## Calculate cell cycle, add to meta data
Using the package scran, get the mouse cell cycle markers and a mapping of m

```r
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
# Convert to matrix for use in cycle
mat <- as.matrix(GetAssayData(experiment.aggregate))

# Convert rownames to ENSEMBL IDs, Using biomaRt
#ensembl<- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#anno_data <- getBM( values=rownames(mat), attributes=c("mgi_symbol","ensembl_gene_id") , filters= "mgi_symbol"  ,mart=ensembl)
# Downloaded from Biomart
#download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/mart_export_June2019.txt", "mart_export_June2019.txt")
anno <- read.delim("mart_export_June2019.txt")

ord <- match(rownames(mat), anno$MGI.symbol) # use anno$mgi_symbol if via biomaRt
rownames(mat) <- anno$Gene.stable.ID[ord] # use anno$ensembl_gene_id if via biomaRt
drop <- which(is.na(rownames(mat)))
mat <- mat[-drop,]
cycles <- scran::cyclone(mat, pairs=mm.pairs)
tmp <- data.frame(cell.cycle = cycles$phases)
rownames(tmp) <- colnames(mat)
experiment.aggregate <- AddMetaData(experiment.aggregate, tmp)
```

### Lets create a fake batch metadata (used in part 3), Here we determine UCD_Adj_VitE is from one batch and UCD_Adj_VitE/UCD_Adj_VitE are from a second battch

Here we build a new metadata variable 'batchid' which can be used to specify treatment groups.

```r
samplename = experiment.aggregate@meta.data$orig.ident
table(samplename)
```

```
## samplename
##  UCD_Adj_VitE UCD_Supp_VitE  UCD_VitE_Def 
##           904          1000           992
```

```r
batchid = rep("Batch1",length(samplename))
batchid[samplename %in% c("UCD_Adj_VitE")] = "Batch2"
names(batchid) = rownames(experiment.aggregate@meta.data)

experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = batchid,
  col.name = "batchid")

table(experiment.aggregate$batchid)
```

```
## 
## Batch1 Batch2 
##   1992    904
```

### Lets spend a little time getting to know the Seurat object.

The Seurat object is the center of each single cell analysis. It stores __all__ information associated with the dataset, including data, annotations, analyses, etc. The R function slotNames can be used to view the slot names within an object.


```r
slotNames(experiment.aggregate)
```

```
##  [1] "assays"       "meta.data"    "active.assay" "active.ident"
##  [5] "graphs"       "neighbors"    "reductions"   "project.name"
##  [9] "misc"         "version"      "commands"     "tools"
```

We can then view the data within a slot with the `@` operator.

```r
head(experiment.aggregate@meta.data)
```

```
##                                  orig.ident nCount_RNA nFeature_RNA
## ACTCTAATGTGGGTATG-UCD_Adj_VitE UCD_Adj_VitE       8885         2818
## AGGCTGGTCAATCACAC-UCD_Adj_VitE UCD_Adj_VitE       5019         2167
## ATGACTAGCACATGACT-UCD_Adj_VitE UCD_Adj_VitE       2208         1286
## AAGCGTCGTCTCTAAGG-UCD_Adj_VitE UCD_Adj_VitE       2795         1474
## ACATCGGGTCCATGCTC-UCD_Adj_VitE UCD_Adj_VitE       5372         2271
## ATACGGTAGTGACCAAG-UCD_Adj_VitE UCD_Adj_VitE        598          367
##                                percent.mito cell.cycle batchid
## ACTCTAATGTGGGTATG-UCD_Adj_VitE     1.969612         G1  Batch2
## AGGCTGGTCAATCACAC-UCD_Adj_VitE     6.216378         G1  Batch2
## ATGACTAGCACATGACT-UCD_Adj_VitE     6.838768         G1  Batch2
## AAGCGTCGTCTCTAAGG-UCD_Adj_VitE     4.221825         G1  Batch2
## ACATCGGGTCCATGCTC-UCD_Adj_VitE     7.557707         G1  Batch2
## ATACGGTAGTGACCAAG-UCD_Adj_VitE    11.371237         G1  Batch2
```

#### Question(s)

1. What slots are empty, what slots have data?
2. What columns are available in meta.data?
3. Look up the help documentation for subset?

## Finally, save the original object, write out a tab-delimited table that could be read into excel, and view the object.

```r
write.table(as.matrix(experiment.data),"raw.datatable.txt",sep="\t",col.names=T,row.names=T)
experiment.aggregate
```

```
## An object of class Seurat 
## 12811 features across 2896 samples within 1 assay 
## Active assay: RNA (12811 features)
```

```r
## Original dataset in Seurat class, with no filtering
save(experiment.aggregate,file="original_seurat_object.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/scRNA_Workshop-PART2.Rmd", "scRNA_Workshop-PART2.Rmd")
```

## Session Information

```r
sessionInfo()
```

```
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Mojave 10.14.5
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] Seurat_3.0.2
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0            Rtsne_0.15                 
##   [3] colorspace_1.4-1            ggridges_0.5.1             
##   [5] dynamicTreeCut_1.63-1       XVector_0.24.0             
##   [7] GenomicRanges_1.36.0        BiocNeighbors_1.2.0        
##   [9] listenv_0.7.0               npsurv_0.4-0               
##  [11] ggrepel_0.8.1               codetools_0.2-16           
##  [13] splines_3.6.0               R.methodsS3_1.7.1          
##  [15] lsei_1.2-0                  knitr_1.23                 
##  [17] scater_1.12.2               jsonlite_1.6               
##  [19] ica_1.0-2                   cluster_2.1.0              
##  [21] png_0.1-7                   R.oo_1.22.0                
##  [23] sctransform_0.2.0           compiler_3.6.0             
##  [25] httr_1.4.0                  dqrng_0.2.1                
##  [27] assertthat_0.2.1            Matrix_1.2-17              
##  [29] lazyeval_0.2.2              limma_3.40.2               
##  [31] BiocSingular_1.0.0          htmltools_0.3.6            
##  [33] tools_3.6.0                 rsvd_1.0.1                 
##  [35] igraph_1.2.4.1              gtable_0.3.0               
##  [37] glue_1.3.1                  GenomeInfoDbData_1.2.1     
##  [39] RANN_2.6.1                  reshape2_1.4.3             
##  [41] dplyr_0.8.1                 Rcpp_1.0.1                 
##  [43] Biobase_2.44.0              gdata_2.18.0               
##  [45] ape_5.3                     nlme_3.1-140               
##  [47] DelayedMatrixStats_1.6.0    gbRd_0.4-11                
##  [49] lmtest_0.9-37               xfun_0.7                   
##  [51] stringr_1.4.0               globals_0.12.4             
##  [53] irlba_2.3.3                 gtools_3.8.1               
##  [55] statmod_1.4.32              future_1.13.0              
##  [57] edgeR_3.26.5                MASS_7.3-51.4              
##  [59] zlibbioc_1.30.0             zoo_1.8-6                  
##  [61] scales_1.0.0                parallel_3.6.0             
##  [63] SummarizedExperiment_1.14.0 RColorBrewer_1.1-2         
##  [65] SingleCellExperiment_1.6.0  yaml_2.2.0                 
##  [67] reticulate_1.12             pbapply_1.4-0              
##  [69] gridExtra_2.3               ggplot2_3.2.0              
##  [71] stringi_1.4.3               S4Vectors_0.22.0           
##  [73] scran_1.12.1                caTools_1.17.1.2           
##  [75] BiocGenerics_0.30.0         BiocParallel_1.18.0        
##  [77] bibtex_0.4.2                GenomeInfoDb_1.20.0        
##  [79] Rdpack_0.11-0               SDMTools_1.1-221.1         
##  [81] rlang_0.3.4                 pkgconfig_2.0.2            
##  [83] bitops_1.0-6                matrixStats_0.54.0         
##  [85] evaluate_0.14               lattice_0.20-38            
##  [87] ROCR_1.0-7                  purrr_0.3.2                
##  [89] htmlwidgets_1.3             cowplot_0.9.4              
##  [91] tidyselect_0.2.5            plyr_1.8.4                 
##  [93] magrittr_1.5                R6_2.4.0                   
##  [95] IRanges_2.18.1              gplots_3.0.1.1             
##  [97] DelayedArray_0.10.0         pillar_1.4.1               
##  [99] fitdistrplus_1.0-14         survival_2.44-1.1          
## [101] RCurl_1.95-4.12             tibble_2.1.3               
## [103] future.apply_1.3.0          tsne_0.1-3                 
## [105] crayon_1.3.4                KernSmooth_2.23-15         
## [107] plotly_4.9.0                rmarkdown_1.13             
## [109] viridis_0.5.1               locfit_1.5-9.1             
## [111] grid_3.6.0                  data.table_1.12.2          
## [113] metap_1.1                   digest_0.6.19              
## [115] tidyr_0.8.3                 R.utils_2.9.0              
## [117] stats4_3.6.0                munsell_0.5.0              
## [119] beeswarm_0.2.3              viridisLite_0.3.0          
## [121] vipor_0.4.5
```
