---
title: "Single Cell RNAseq Part 3"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---
## Load libraries

```r
library(Seurat)
```

## Load the Seurat object

```r
load(file="pre_sample_corrected.RData")
experiment.aggregate
```

```
## An object of class Seurat 
## 12811 features across 2681 samples within 1 assay 
## Active assay: RNA (12811 features)
```

## Exploring Batch effects 3 ways, none, Seurat [vars.to.regress] and COMBAT

First lets view the data without any corrections

## PCA in prep for tSNE

ScaleData - Scales and centers genes in the dataset. 

```r
# ?ScaleData

experiment.aggregate.noc <- ScaleData(object = experiment.aggregate)
```

```
## Centering and scaling data matrix
```

Run PCA

```r
experiment.aggregate.noc <- RunPCA(object = experiment.aggregate.noc)
```

```
## PC_ 1 
## Positive:  Txn1, Sncg, Fez1, S100a10, Atp6v0b, Lxn, Dctn3, Tppp3, Sh3bgrl3, Rabac1 
## 	   Cisd1, Ppia, Fxyd2, Atp6v1f, Stmn3, Ndufa11, Atp5g1, Bex2, Atpif1, Uchl1 
## 	   Ndufa4, Psmb6, Ubb, Hagh, Anxa2, Gabarapl2, Rgs10, Nme1, Prdx2, Psmb3 
## Negative:  Mt1, Malat1, Adcyap1, Ptn, Apoe, Zeb2, Mt2, Timp3, Fabp7, Gpm6b 
## 	   Gal, Kit, Qk, Atp2b4, Plp1, Ifitm3, 6330403K07Rik, Sparc, Id3, Gap43 
## 	   Selenop, Gpx3, Zfp36l1, Rgcc, Scg2, Cbfb, Zfp36, Igfbp7, Marcksl1, Phlda1 
## PC_ 2 
## Positive:  Nefh, Cntn1, Thy1, S100b, Sv2b, Cplx1, Slc17a7, Vamp1, Nefm, Lynx1 
## 	   Endod1, Atp1b1, Scn1a, Vsnl1, Nat8l, Ntrk3, Sh3gl2, Fam19a2, Eno2, Scn1b 
## 	   Spock1, Scn8a, Glrb, Syt2, Scn4b, Lrrn1, Lgi3, Atp2b2, Cpne6, Snap25 
## Negative:  Malat1, Cd24a, Tmem233, Cd9, Dusp26, Mal2, Tmem158, Carhsp1, Fxyd2, Ctxn3 
## 	   Prkca, Ubb, Arpc1b, Crip2, Gna14, S100a6, Cd44, Tmem45b, Klf5, Tceal9 
## 	   Hs6st2, Bex3, Cd82, Emp3, Ift122, Tubb2b, Pfn1, Fam89a, Dynll1, Acpp 
## PC_ 3 
## Positive:  P2ry1, Fam19a4, Gm7271, Rarres1, Th, Zfp521, Wfdc2, Tox3, Gfra2, D130079A08Rik 
## 	   Iqsec2, Pou4f2, Rgs5, Kcnd3, Id4, Rasgrp1, Slc17a8, Casz1, Cdkn1a, Piezo2 
## 	   Dpp10, Gm11549, Fxyd6, Spink2, Rgs10, Zfhx3, C1ql4, Cd34, Gabra1, Cckar 
## Negative:  Calca, Basp1, Gap43, Ppp3ca, Map1b, Scg2, Cystm1, Tmem233, Map7d2, Calm1 
## 	   Ift122, Tubb3, Ncdn, Resp18, Prkca, Etv1, Nmb, Skp1a, Crip2, Camk2a 
## 	   Epb41l3, Tspan8, Ntrk1, Gna14, Deptor, Adk, Jak1, Tmem255a, Etl4, Camk2g 
## PC_ 4 
## Positive:  Id3, Timp3, Selenop, Pvalb, Ifitm3, Sparc, Igfbp7, Adk, Tm4sf1, Sgk1 
## 	   Ly6c1, Id1, Etv1, Nsg1, Mt2, Slc17a7, Zfp36l1, Spp1, Cldn5, Itm2a 
## 	   Ier2, Aldoc, Shox2, Cxcl12, Ptn, Stxbp6, Qk, Sparcl1, Jak1, Slit2 
## Negative:  Gap43, Calca, Stmn1, Tac1, Ppp3ca, 6330403K07Rik, Arhgdig, Alcam, Adcyap1, Prune2 
## 	   Kit, Ngfr, Ywhag, Gal, Atp1a1, Fxyd6, Smpd3, Ntrk1, Tmem100, Mt3 
## 	   Atp2b4, Cd24a, Cnih2, Tppp3, Gpx3, S100a11, Scn7a, Snap25, Cbfb, Gnb1 
## PC_ 5 
## Positive:  Cpne3, Klf5, Acpp, Fxyd2, Jak1, Nppb, Osmr, Rgs4, Htr1f, Nbl1 
## 	   Gm525, Etv1, Sst, Zfhx3, Adk, Tspan8, Cysltr2, Parm1, Npy2r, Tmem233 
## 	   Nts, Cd24a, Prkca, Prune2, Socs2, Il31ra, Dgkz, Gnb1, Phf24, Plxnc1 
## Negative:  Mt1, Ptn, B2m, Prdx1, Dbi, Mt3, Ifitm3, Fxyd7, Id3, S100a16 
## 	   Calca, Mt2, Sparc, Pcp4l1, Selenop, Ifitm2, Rgcc, Igfbp7, Abcg2, Tm4sf1 
## 	   Selenom, S100a13, Apoe, Hspb1, Timp3, Cryab, Ubb, Ly6c1, Phlda1, Gap43
```

```r
DimPlot(object = experiment.aggregate.noc, reduction = "pca")
```

<img src="scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

PCA Elbow plot to determine how many principal components to use in downstream analyses.  Components after the "elbow" in the plot generally explain little additional variability in the data.


```r
ElbowPlot(experiment.aggregate.noc)
```

![](scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

We use 10 components in downstream analyses. Using more components more closely approximates the full data set but increases run time.

TSNE Plot

```r
pcs.use <- 10
experiment.aggregate.noc <- RunTSNE(object = experiment.aggregate.noc, dims = 1:pcs.use)
DimPlot(object = experiment.aggregate.noc)
```

<img src="scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

## Correct for sample to sample differences (seurat)

Use vars.to.regress to correct for the sample to sample differences and percent mitochondria

```r
experiment.aggregate.regress <- ScaleData(object = experiment.aggregate, 
                                          vars.to.regress = c("orig.ident", "percent.mito"), model.use = "linear")
```

```
## Regressing out orig.ident, percent.mito
```

```
## Centering and scaling data matrix
```

```r
experiment.aggregate.regress <- RunPCA(object =experiment.aggregate.regress)
```

```
## PC_ 1 
## Positive:  Txn1, Sncg, Fxyd2, Lxn, Fez1, Rgs10, Synpr, S100a10, Atp6v0b, Ppp1r1a 
## 	   Dctn3, Ubb, Tmem45b, Rabac1, Sh3bgrl3, Tppp3, Atpif1, Bex2, Cisd1, Atp6v1f 
## 	   Bex3, Hagh, Tmsb4x, Psmb6, Ndufa11, Pfn1, Aldoa, Ndufs5, Anxa2, Ppia 
## Negative:  Ptn, S100b, Cbfb, Mt1, Sv2b, Timp3, Nfia, Lynx1, Map2, Ngfr 
## 	   Adcyap1, Fxyd7, Gap43, Nefh, Thy1, Enah, Syt2, Scg2, Igfbp7, Tmem229b 
## 	   Faim2, Nptn, Kit, Nfib, Slc17a7, Ryr2, Epb41l3, Cntnap2, Vat1l, Zeb2 
## PC_ 2 
## Positive:  Cntn1, Nefh, Cplx1, Thy1, Slc17a7, S100b, Sv2b, Ntrk3, Atp1b1, Scn1a 
## 	   Nefm, Vamp1, Lrrn1, Atp2b2, Hopx, Endod1, Tagln3, Scn1b, Snap25, Vsnl1 
## 	   Nat8l, Nefl, Lynx1, Glrb, Scn4b, Eno2, Fam19a2, Sh3gl2, Scn8a, Spock1 
## Negative:  Malat1, Tmem233, Cd9, Cd24a, Prkca, Mal2, Dusp26, Carhsp1, Gna14, Crip2 
## 	   Gadd45g, Id3, Osmr, Ift122, Tmem158, Cd44, Calca, Camk2a, Hs6st2, Cd82 
## 	   Emp3, Gm525, Ctxn3, Crip1, Nppb, S100a6, Socs2, Sst, Tac1, Mt1 
## PC_ 3 
## Positive:  P2ry1, Fam19a4, Gm7271, Rarres1, Th, Zfp521, Wfdc2, Cdkn1a, Tox3, Gfra2 
## 	   D130079A08Rik, Rgs5, Cd81, Kcnd3, Pou4f2, Cd34, Iqsec2, Sorbs2, Slc17a8, Rasgrp1 
## 	   Casz1, Id4, Dpp10, Piezo2, Igfbp7, Gm11549, Zfhx3, Spink2, Ifitm3, Gabra1 
## Negative:  Calca, Ppp3ca, Basp1, Map1b, Gap43, Cystm1, Tubb3, Scg2, Calm1, Ncdn 
## 	   Map7d2, Epb41l3, Ift122, 6330403K07Rik, Skp1a, Tmem233, Nmb, Dusp26, Tmem255a, Crip2 
## 	   Resp18, Ntrk1, Ywhag, Prkca, Tnfrsf21, Camk2a, Fxyd7, Deptor, Mt3, Gna14 
## PC_ 4 
## Positive:  Adk, Pvalb, Etv1, Nsg1, Jak1, Tmem233, Tspan8, Nppb, Slc17a7, Shox2 
## 	   Sst, Gm525, Htr1f, Spp1, Slit2, Cbln2, Nts, Stxbp6, Aldoc, Osmr 
## 	   Cmtm8, Runx3, Fam19a2, Cysltr2, Hapln4, Rasgrp2, Aprt, Carhsp1, Nxph1, Atp1a3 
## Negative:  Gap43, Calca, Arhgdig, Stmn1, Alcam, Tac1, Ngfr, 6330403K07Rik, Ppp3ca, Kit 
## 	   Fxyd6, Adcyap1, Smpd3, Atp1a1, Ntrk1, Gm7271, Tagln3, Gal, Tmem100, Prune2 
## 	   Chl1, Atp2b4, Tppp3, Dclk1, Cnih2, Mgll, S100a11, Scn7a, Mt3, Bcl11b 
## PC_ 5 
## Positive:  Fxyd2, Rgs4, Cpne3, Acpp, Klf5, Zfhx3, Nbl1, Prune2, Cd24a, Phf24 
## 	   Gnb1, Osmr, Jak1, Dgkz, Prkca, Parm1, Tmem233, Tspan8, Nppb, Etv1 
## 	   Kif5b, Socs2, Synpr, Dpp10, Plxnc1, Ywhag, Htr1f, Gm525, Adk, Cysltr2 
## Negative:  Prdx1, Mt1, Dbi, Ptn, B2m, Mt3, Ubb, Id3, Sparc, Mt2 
## 	   Ifitm3, Fxyd7, Selenop, Rgcc, Tecr, Cryab, Apoe, Dad1, Uqcrb, Hspa1a 
## 	   Phlda1, Spcs1, Fxyd1, Selenom, Psmb2, Timp3, S100a16, Qk, Mgst3, Gpm6b
```

```r
DimPlot(object = experiment.aggregate.regress, reduction.use = "pca")
```

<img src="scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

Corrected TSNE Plot

```r
experiment.aggregate.regress <- RunTSNE(object = experiment.aggregate.regress, dims.use = 1:pcs.use)
DimPlot(object = experiment.aggregate.regress, reduction = "tsne")
```

<img src="scRNA_Workshop-PART3_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

## COMBAT corrected, https://academic.oup.com/biostatistics/article-lookup/doi/10.1093/biostatistics/kxj037

## NOT RUN

```r
# Install sva package if not installed
# For R 3.4 and earlier use
# source("https://bioconductor.org/biocLite.R")
# biocLite("sva")
#
# For R 3.5 and later use
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("sva")
library(sva)
?ComBat
m = as.data.frame(as.matrix(GetAssayData(experiment.aggregate)))
com = ComBat(dat=m, batch=as.numeric(as.factor(experiment.aggregate$orig.ident)), prior.plots=FALSE, par.prior=TRUE)
```



```r
experiment.aggregate.combat <- experiment.aggregate
SetAssayData(experiment.aggregate.combat, new.data = as.matrix(com))
experiment.aggregate.combat = ScaleData(experiment.aggregate.combat)
```

Principal components on ComBat adjusted data

```r
experiment.aggregate.combat <- RunPCA(object = experiment.aggregate.combat)

DimPlot(object = experiment.aggregate.combat, reduction = "pca")
```

TSNE plot on ComBat adjusted data

```r
experiment.aggregate.combat <- RunTSNE(object = experiment.aggregate.combat, dims.use = 1:pcs.use)
DimPlot(object = experiment.aggregate.combat, reduction = "tsne")
```

#### Question(s)

1. Try a couple of PCA cutoffs (low and high) and compare the TSNE plots from the different methods.  Do they look meaningfully different?

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/scRNA_Workshop-PART4.Rmd", "scRNA_Workshop-PART4.Rmd")
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
## [67] igraph_1.2.4.1      labeling_0.3        bitops_1.0-6       
## [70] rmarkdown_1.13      npsurv_0.4-0        gtable_0.3.0       
## [73] codetools_0.2-16    reshape2_1.4.3      R6_2.4.0           
## [76] gridExtra_2.3       zoo_1.8-6           knitr_1.23         
## [79] dplyr_0.8.1         future.apply_1.3.0  KernSmooth_2.23-15 
## [82] metap_1.1           ape_5.3             stringi_1.4.3      
## [85] parallel_3.6.0      Rcpp_1.0.1          sctransform_0.2.0  
## [88] png_0.1-7           tidyselect_0.2.5    xfun_0.7           
## [91] lmtest_0.9-37
```
