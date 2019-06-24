---
title: "Single Cell RNAseq Part 1"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---
[Seurat](http://satijalab.org/seurat/) is a popular R package that is designed for QC, analysis, and exploration of single cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single cell transcriptomic measurements, and to integrate diverse types of single cell data. Further, the authors provide several [tutorials](http://satijalab.org/seurat/get_started.html), on their website.

Dowload and expand the expression_tables.tar.gz file to extract the single cell matrix files for the three samples. These are isolated mouse cells ran on the 10X genomics platform for single cell RNA sequencing, sequenced with UC Davis on 1 HiSeq 4000.

* UCD_VitE_Def
* UCD_Supp_VitE
* UCD_Adj_VitE

We start with loading needed libraries for R, at this time all we need is the package [Seurat](http://satijalab.org/seurat/).

```r
library(Seurat)
```

## Load the Cell Ranger Matrix Data and create the base Seurat object.
Cell Ranger provides a function `cellranger aggr` that will combine multiple samples into a single matrix file. However, when processing data in R and Seurat this is unnecessary and we can aggregate them in R.

Seurat provides a function `Read10X` to read in 10X data folder. First we read in data from each individual sample folder. First, we initialize the Seurat object (`CreateSeuratObject`) with the raw (non-normalized data). Keep all genes expressed in >= 10 cells. Keep all cells with at least 200 detected genes. Also extracting sample names, calculating and adding in the metadata mitochondrial percentage of each cell. Adding in the metadata batchid. Finally, saving the raw Seurat object.

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
### Calc mitocondrial content
Calculate percent mitochondrial genes per cell. In mouse these genes can be identified as those that begin with 'mt', in human data they begin with MT.

```r
mito.genes <- grep("^mt-", rownames(experiment.aggregate), value = T)
percent.mito <- Matrix::colSums(GetAssayData(experiment.aggregate, slot = "counts")[mito.genes, ]) / Matrix::colSums(GetAssayData(experiment.aggregate, slot = "counts"))

# Add to @meta.data
experiment.aggregate$percent.mito <- percent.mito
```
### Lets fix the sample names, reassign names with more meaningful factors

The original samples names (the names above in ids) can be found in the metadata slot, column orig.ident. Here we build a new metadata variable 'batchid' which can be used to specify treatment groups.

```r
samplename = experiment.aggregate$orig.ident
table(samplename)
```

```
## samplename
##  UCD_Adj_VitE UCD_Supp_VitE  UCD_VitE_Def 
##           904          1000           992
```

```r
batchid = experiment.aggregate$orig.ident
names(batchid) = rownames(experiment.aggregate@meta.data)

experiment.aggregate$batchid <- batchid
table(experiment.aggregate$batchid)
```

```
## 
##  UCD_Adj_VitE UCD_Supp_VitE  UCD_VitE_Def 
##           904          1000           992
```
### Lets spend a little time getting to know the Seurat object.

The Seurat object is the center of each single cell analysis. It stores __all__ information associated with the dataset, including data, annotations, analyes, etc. The R function slotNames can be used to view the slot names within an object.


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
##                                percent.mito      batchid
## ACTCTAATGTGGGTATG-UCD_Adj_VitE   0.01969612 UCD_Adj_VitE
## AGGCTGGTCAATCACAC-UCD_Adj_VitE   0.06216378 UCD_Adj_VitE
## ATGACTAGCACATGACT-UCD_Adj_VitE   0.06838768 UCD_Adj_VitE
## AAGCGTCGTCTCTAAGG-UCD_Adj_VitE   0.04221825 UCD_Adj_VitE
## ACATCGGGTCCATGCTC-UCD_Adj_VitE   0.07557707 UCD_Adj_VitE
## ATACGGTAGTGACCAAG-UCD_Adj_VitE   0.11371237 UCD_Adj_VitE
```

#### Question(s)

1. What slots are empty, what slots have data?
2. What columns are available in meta.data?
3. Look up the help documentation for subset?

## Finally, save the original object, write out a tab-delimited table that could be read into excel, and view the object.

```r
# write.table(as.matrix(experiment.data),"raw.datatable.txt",sep="\t",col.names=T,row.names=T)
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
##  [1] httr_1.4.0          tidyr_0.8.3         viridisLite_0.3.0  
##  [4] jsonlite_1.6        splines_3.6.0       lsei_1.2-0         
##  [7] R.utils_2.9.0       gtools_3.8.1        Rdpack_0.11-0      
## [10] assertthat_0.2.1    yaml_2.2.0          ggrepel_0.8.1      
## [13] globals_0.12.4      pillar_1.4.1        lattice_0.20-38    
## [16] reticulate_1.12     glue_1.3.1          digest_0.6.19      
## [19] RColorBrewer_1.1-2  SDMTools_1.1-221.1  colorspace_1.4-1   
## [22] cowplot_0.9.4       htmltools_0.3.6     Matrix_1.2-17      
## [25] R.oo_1.22.0         plyr_1.8.4          pkgconfig_2.0.2    
## [28] bibtex_0.4.2        tsne_0.1-3          listenv_0.7.0      
## [31] purrr_0.3.2         scales_1.0.0        RANN_2.6.1         
## [34] gdata_2.18.0        Rtsne_0.15          tibble_2.1.3       
## [37] ggplot2_3.2.0       ROCR_1.0-7          pbapply_1.4-0      
## [40] lazyeval_0.2.2      survival_2.44-1.1   magrittr_1.5       
## [43] crayon_1.3.4        evaluate_0.14       R.methodsS3_1.7.1  
## [46] future_1.13.0       nlme_3.1-140        MASS_7.3-51.4      
## [49] gplots_3.0.1.1      ica_1.0-2           tools_3.6.0        
## [52] fitdistrplus_1.0-14 data.table_1.12.2   gbRd_0.4-11        
## [55] stringr_1.4.0       plotly_4.9.0        munsell_0.5.0      
## [58] cluster_2.1.0       irlba_2.3.3         compiler_3.6.0     
## [61] rsvd_1.0.1          caTools_1.17.1.2    rlang_0.3.4        
## [64] grid_3.6.0          ggridges_0.5.1      htmlwidgets_1.3    
## [67] igraph_1.2.4.1      bitops_1.0-6        rmarkdown_1.13     
## [70] npsurv_0.4-0        gtable_0.3.0        codetools_0.2-16   
## [73] reshape2_1.4.3      R6_2.4.0            gridExtra_2.3      
## [76] zoo_1.8-6           knitr_1.23          dplyr_0.8.1        
## [79] future.apply_1.3.0  KernSmooth_2.23-15  metap_1.1          
## [82] ape_5.3             stringi_1.4.3       parallel_3.6.0     
## [85] Rcpp_1.0.1          sctransform_0.2.0   png_0.1-7          
## [88] tidyselect_0.2.5    xfun_0.7            lmtest_0.9-37
```
