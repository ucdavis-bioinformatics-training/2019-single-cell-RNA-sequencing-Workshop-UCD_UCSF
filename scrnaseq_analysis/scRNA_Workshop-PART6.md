---
title: "Single Cell RNAseq Part 6"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Reminder of samples

* UCD_Adj_VitE
* UCD_Supp_VitE
* UCD_VitE_Def

## Load libraries

```r
library(Seurat)
library(ggplot2)
```

## Load the Seurat object

```r
load("clusters_seurat_object.RData")
experiment.merged
```

```
## An object of class Seurat 
## 12811 features across 2681 samples within 1 assay 
## Active assay: RNA (12811 features)
##  2 dimensional reductions calculated: pca, tsne
```


#0. Setup
Load the final Seurat object, load libraries (also see additional required packages for each example)


#1. DE With Single Cell Data Using Limma
For differential expression using models more complex than those allowed by FindAllMarkers(), data from Seurat may be used in limma (https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf)

We illustrate by comparing sample 1 to sample 2 within cluster 0:

```r
library(limma)
cluster0 <- subset(experiment.merged, idents = '0')
expr <- as.matrix(GetAssayData(cluster0))

# Filter out genes that are 0 for every cell in this cluster
bad <- which(rowSums(expr) == 0)
expr <- expr[-bad,]

mm <- model.matrix(~0 + orig.ident, data = cluster0@meta.data)
fit <- lmFit(expr, mm)  
head(coef(fit)) # means in each sample for each gene
```

```
##        orig.identUCD_Adj_VitE orig.identUCD_Supp_VitE orig.identUCD_VitE_Def
## Xkr4               0.00000000             0.000000000            0.007299386
## Sox17              0.00000000             0.009059615            0.000000000
## Mrpl15             0.09950052             0.047111832            0.090124656
## Lypla1             0.16825245             0.174572188            0.274566585
## Tcea1              0.23337111             0.217451945            0.268055981
## Rgs20              0.05671722             0.098139805            0.075264334
```

```r
contr <- makeContrasts(orig.identUCD_Supp_VitE - orig.identUCD_Adj_VitE, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrasts = contr)
tmp <- eBayes(tmp)
topTable(tmp, sort.by = "P", n = 20) # top 20 DE genes
```

```
##              logFC   AveExpr         t      P.Value    adj.P.Val         B
## Rpl21   -0.6940744 2.0960659 -7.491480 1.695939e-13 2.114836e-09 20.010067
## Rpl23a  -0.5758956 2.1549568 -6.358577 3.308626e-10 2.062928e-06 12.817327
## Tmsb10  -0.4870359 3.2805433 -6.263147 5.964289e-10 2.479156e-06 12.259827
## Rpl39   -0.5102379 2.2270006 -6.138103 1.275859e-09 3.977489e-06 11.540934
## Rpl17   -0.5555388 1.9858858 -5.966323 3.548324e-09 8.849521e-06 10.574973
## Rpl32   -0.5305218 2.4879312 -5.922604 4.584875e-09 9.528899e-06 10.333143
## Rps8    -0.5176911 2.4850790 -5.631597 2.420869e-08 4.120323e-05  8.765247
## S100a10 -0.5327206 1.2593306 -5.615855 2.643351e-08 4.120323e-05  8.682515
## Rps11   -0.4909334 2.2978636 -5.559438 3.616049e-08 5.004915e-05  8.387774
## Tmsb4x  -0.3659831 3.3269827 -5.529852 4.257135e-08 5.004915e-05  8.234309
## Pcp4    -0.6658739 1.3585640 -5.523235 4.414921e-08 5.004915e-05  8.200095
## Hspa8   -0.5093576 2.1252230 -5.445758 6.741872e-08 7.005928e-05  7.802273
## Rps24   -0.4879797 1.5510296 -5.303590 1.446206e-07 1.380947e-04  7.085914
## Rpl23   -0.4631348 2.0116972 -5.290467 1.550382e-07 1.380947e-04  7.020680
## Sncg    -0.4462357 3.2198350 -5.273893 1.692372e-07 1.406925e-04  6.938511
## Pfn1    -0.4909951 1.3084847 -5.246140 1.958794e-07 1.526635e-04  6.801457
## Rpl24   -0.4670874 2.0841996 -5.193872 2.574962e-07 1.888810e-04  6.545176
## H3f3b   -0.4806824 2.0994168 -5.167956 2.946317e-07 2.041143e-04  6.418997
## Car8    -0.1846817 0.1153888 -5.139171 3.419529e-07 2.170852e-04  6.279540
## Fau     -0.3994305 2.7479232 -5.135678 3.481720e-07 2.170852e-04  6.262668
```
* logFC: log2 fold change (UCD_Supp_VitE/UCD_Adj_VitE)
* AveExpr: Average expression, in log2 counts per million, across all cells included in analysis (i.e. those in cluster 0)
* t: t-statistic, i.e. logFC divided by its standard error
* P.Value: Raw p-value from test that logFC differs from 0
* adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value

The limma vignette linked above gives more detail on model specification.

# 2. Gene Ontology (GO) Enrichment of Genes Expressed in a Cluster

```r
library(topGO)
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport, clusterMap, parApply, parCapply, parLapply, parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply, union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: graph
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with 'browseVignettes()'. To cite Bioconductor, see 'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: GO.db
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## 
```

```
## Loading required package: SparseM
```

```
## 
## Attaching package: 'SparseM'
```

```
## The following object is masked from 'package:base':
## 
##     backsolve
```

```
## 
## groupGOTerms: 	GOBPTerm, GOMFTerm, GOCCTerm environments built.
```

```
## 
## Attaching package: 'topGO'
```

```
## The following object is masked from 'package:IRanges':
## 
##     members
```

```r
# install org.Mm.eg.db from Bioconductor if not already installed (for mouse only)
cluster0 <- subset(experiment.merged, idents = '0')
expr <- as.matrix(GetAssayData(cluster0))
# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.75)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
	GOdata <- new("topGOdata",
		ontology = "BP", # use biological process ontology
		allGenes = geneList,
		geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
```

```
## 
## Building most specific GOs .....
```

```
## Loading required package: org.Mm.eg.db
```

```
## 
```

```
## 	( 10570 GO terms found. )
```

```
## 
## Build GO DAG topology ..........
```

```
## 	( 14604 GO terms and 34657 relations. )
```

```
## 
## Annotating nodes ...............
```

```
## 	( 11727 genes annotated to the GO terms. )
```

```r
# Test for enrichment using Fisher's Exact Test
	resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
```

```
## 
## 			 -- Elim Algorithm -- 
## 
## 		 the algorithm is scoring 3146 nontrivial nodes
## 		 parameters: 
## 			 test statistic: fisher
## 			 cutOff: 0.01
```

```
## 
## 	 Level 19:	1 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 18:	1 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 17:	2 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 16:	6 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 15:	17 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 14:	34 nodes to be scored	(31 eliminated genes)
```

```
## 
## 	 Level 13:	82 nodes to be scored	(38 eliminated genes)
```

```
## 
## 	 Level 12:	125 nodes to be scored	(438 eliminated genes)
```

```
## 
## 	 Level 11:	201 nodes to be scored	(521 eliminated genes)
```

```
## 
## 	 Level 10:	290 nodes to be scored	(750 eliminated genes)
```

```
## 
## 	 Level 9:	358 nodes to be scored	(999 eliminated genes)
```

```
## 
## 	 Level 8:	409 nodes to be scored	(1361 eliminated genes)
```

```
## 
## 	 Level 7:	496 nodes to be scored	(1659 eliminated genes)
```

```
## 
## 	 Level 6:	482 nodes to be scored	(2334 eliminated genes)
```

```
## 
## 	 Level 5:	340 nodes to be scored	(2567 eliminated genes)
```

```
## 
## 	 Level 4:	194 nodes to be scored	(2588 eliminated genes)
```

```
## 
## 	 Level 3:	88 nodes to be scored	(2658 eliminated genes)
```

```
## 
## 	 Level 2:	19 nodes to be scored	(3206 eliminated genes)
```

```
## 
## 	 Level 1:	1 nodes to be scored	(3206 eliminated genes)
```

```r
	GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
```

```
##         GO.ID                                                            Term Annotated Significant Expected  Fisher
## 1  GO:0006412                                                     translation       517          52     6.83 1.0e-22
## 2  GO:0002181                                         cytoplasmic translation        77          18     1.02 4.7e-18
## 3  GO:0000028                                ribosomal small subunit assembly        18           7     0.24 1.7e-09
## 4  GO:0000027                                ribosomal large subunit assembly        29           7     0.38 7.5e-08
## 5  GO:0097214          positive regulation of lysosomal membrane permeability         2           2     0.03 0.00017
## 6  GO:0006880                          intracellular sequestering of iron ion         3           2     0.04 0.00052
## 7  GO:2000582 positive regulation of ATP-dependent microtubule motor activ...         3           2     0.04 0.00052
## 8  GO:0000462 maturation of SSU-rRNA from tricistronic rRNA transcript (SS...        30           4     0.40 0.00062
## 9  GO:0007409                                                    axonogenesis       353          15     4.67 0.00065
## 10 GO:1990090               cellular response to nerve growth factor stimulus        32           4     0.42 0.00079
## 11 GO:0016198                                   axon choice point recognition         4           2     0.05 0.00102
## 12 GO:0050848                        regulation of calcium-mediated signaling        60           5     0.79 0.00115
## 13 GO:0061844 antimicrobial humoral immune response mediated by antimicrob...        17           3     0.22 0.00134
## 14 GO:0050774                   negative regulation of dendrite morphogenesis        18           3     0.24 0.00160
## 15 GO:0002227                                innate immune response in mucosa         5           2     0.07 0.00169
## 16 GO:0071635 negative regulation of transforming growth factor beta produ...         5           2     0.07 0.00169
## 17 GO:0047497                       mitochondrion transport along microtubule        19           3     0.25 0.00188
## 18 GO:0048588                                       developmental cell growth       213           9     2.82 0.00207
## 19 GO:0033138          positive regulation of peptidyl-serine phosphorylation        71           5     0.94 0.00244
## 20 GO:0002679                  respiratory burst involved in defense response         6           2     0.08 0.00251
```
* Annotated: number of genes (out of all.genes) that are annotated with that GO term
* Significant: number of genes that are annotated with that GO term and meet our criteria for "expressed"
* Expected: Under random chance, number of genes that would be expected to be annotated with that GO term and meeting our criteria for "expressed"
* Fisher: (Raw) p-value from Fisher's Exact Test

#3. Weighted Gene Co-Expression Network Analysis (WGCNA)
WGCNA identifies groups of genes ("modules") with correlated expression.
WARNING: TAKES A LONG TIME TO RUN

```r
library(WGCNA)
```

```
## Loading required package: dynamicTreeCut
```

```
## Loading required package: fastcluster
```

```
## 
## Attaching package: 'fastcluster'
```

```
## The following object is masked from 'package:stats':
## 
##     hclust
```

```
## 
## Attaching package: 'WGCNA'
```

```
## The following object is masked from 'package:IRanges':
## 
##     cor
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     cor
```

```
## The following object is masked from 'package:stats':
## 
##     cor
```

```r
options(stringsAsFactors = F)
datExpr <- t(as.matrix(GetAssayData(experiment.merged)))[,VariableFeatures(experiment.merged)]  # only use variable genes in analysis

net <- blockwiseModules(datExpr, power = 10,
  corType = "bicor", # use robust correlation
	networkType = "signed", minModuleSize = 10,
	reassignThreshold = 0, mergeCutHeight = 0.15,
	numericLabels = F, pamRespectsDendro = FALSE,
	saveTOMs = TRUE,
	saveTOMFileBase = "TOM",
	verbose = 3)
```

```
##  Calculating module eigengenes block-wise from all genes
##    Flagging genes and samples with too many missing values...
##     ..step 1
##  ..Working on block 1 .
##     TOM calculation: adjacency..
##     ..will not use multithreading.
## alpha: 1.000000
##      Fraction of slow calculations: 0.000000
##     ..connectivity..
##     ..matrix multiplication (system BLAS)..
##     ..normalization..
##     ..done.
##    ..saving TOM for block 1 into file TOM-block.1.RData
##  ....clustering..
##  ....detecting modules..
##  ....calculating module eigengenes..
##  ....checking kME in modules..
##      ..removing 67 genes from module 1 because their KME is too low.
##      ..removing 43 genes from module 3 because their KME is too low.
##      ..removing 2 genes from module 12 because their KME is too low.
##  ..merging modules that are too close..
##      mergeCloseModules: Merging modules whose distance is less than 0.15
## alpha: 1.000000
##        Calculating new MEs...
## alpha: 1.000000
```

```r
table(net$colors)
```

```
## 
##     black      blue     brown     green      grey       red turquoise    yellow 
##        11        80        21        12      1536        11       312        17
```

```r
# Convert labels to colors for plotting
mergedColors = net$colors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```

![](scRNA_Workshop-PART6_files/figure-html/unnamed-chunk-6-1.png)<!-- -->
Genes in grey module are unclustered.

What genes are in the "blue" module?

```r
colnames(datExpr)[net$colors == "blue"]
```

```
##  [1] "Lxn"           "Txn1"          "Grik1"         "Fez1"          "Tmem45b"       "Synpr"         "Tceal9"        "Ppp1r1a"       "Rgs10"         "Nrn1"          "Fxyd2"         "Ostf1"         "Lix1"          "Sncb"          "Paqr5"         "Bex3"          "Anxa5"         "Gfra2"         "Scg3"          "Ppm1j"         "Kcnab1"        "Kcnip4"        "Cadm1"         "Isl2"          "Pla2g7"        "Tppp3"         "Rgs4"         
## [28] "Tmsb4x"        "Unc119"        "Pmm1"          "Ccdc68"        "Rnf7"          "Prr13"         "Rsu1"          "Pmp22"         "Acpp"          "Kcnip2"        "Cdk15"         "Mrps6"         "Ebp"           "Hexb"          "Cdh11"         "Dapk2"         "Ano3"          "Pde6d"         "Snx7"          "Dtnbp1"        "Tubb2b"        "Nr2c2ap"       "Phf24"         "Rcan2"         "Fam241b"       "Pmvk"          "Slc25a4"      
## [55] "Zfhx3"         "Dgkz"          "Ndufv1"        "Ptrh1"         "1700037H04Rik" "Kif5b"         "Sae1"          "Sri"           "Cpne3"         "Dgcr6"         "Cisd3"         "Syt7"          "Lhfpl3"        "Dda1"          "Ppp1ca"        "Glrx3"         "Stoml1"        "Plagl1"        "Lbh"           "Degs1"         "AI413582"      "Car10"         "Tlx2"          "Parm1"         "March11"       "Cpe"
```

Each cluster is represented by a summary "eigengene".
Plot eigengenes for each non-grey module by clusters from Seurat:

```r
f <- function(module){
  eigengene <- unlist(net$MEs[paste0("ME", module)])
  means <- tapply(eigengene, Idents(experiment.merged), mean, na.rm = T)
  return(means)
}
modules <- c("blue", "brown", "green", "turquoise", "yellow")
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = "Seurat Cluster",
        ylab = "WGCNA Module Eigengene")
axis(1, at = 1:19, labels = 0:18, cex.axis = 0.8)
matpoints(plotdat, col = modules, pch = 21)
```

![](scRNA_Workshop-PART6_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/scRNA_Workshop-PART7.Rmd", "scRNA_Workshop-PART7.Rmd")
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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] WGCNA_1.68            fastcluster_1.1.25    dynamicTreeCut_1.63-1 org.Mm.eg.db_3.8.2    topGO_2.36.0          SparseM_1.77          GO.db_3.8.2           AnnotationDbi_1.46.0  IRanges_2.18.1        S4Vectors_0.22.0      Biobase_2.44.0        graph_1.62.0          BiocGenerics_0.30.0   limma_3.40.2          ggplot2_3.2.0         Seurat_3.0.2         
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.15            colorspace_1.4-1      ggridges_0.5.1        htmlTable_1.13.1      base64enc_0.1-3       rstudioapi_0.10       listenv_0.7.0         npsurv_0.4-0          ggrepel_0.8.1         bit64_0.9-7           mvtnorm_1.0-11        codetools_0.2-16      splines_3.6.0         R.methodsS3_1.7.1     doParallel_1.0.14     robustbase_0.93-5     impute_1.58.0         lsei_1.2-0            knitr_1.23            Formula_1.2-3        
##  [21] jsonlite_1.6          ica_1.0-2             cluster_2.1.0         png_0.1-7             R.oo_1.22.0           sctransform_0.2.0     rrcov_1.4-7           compiler_3.6.0        httr_1.4.0            backports_1.1.4       assertthat_0.2.1      Matrix_1.2-17         lazyeval_0.2.2        acepack_1.4.1         htmltools_0.3.6       tools_3.6.0           rsvd_1.0.1            igraph_1.2.4.1        gtable_0.3.0          glue_1.3.1           
##  [41] RANN_2.6.1            reshape2_1.4.3        dplyr_0.8.1           Rcpp_1.0.1            preprocessCore_1.46.0 gdata_2.18.0          ape_5.3               nlme_3.1-140          iterators_1.0.10      gbRd_0.4-11           lmtest_0.9-37         xfun_0.7              stringr_1.4.0         globals_0.12.4        irlba_2.3.3           gtools_3.8.1          future_1.13.0         DEoptimR_1.0-8        MASS_7.3-51.4         zoo_1.8-6            
##  [61] scales_1.0.0          RColorBrewer_1.1-2    yaml_2.2.0            memoise_1.1.0         reticulate_1.12       pbapply_1.4-0         gridExtra_2.3         rpart_4.1-15          latticeExtra_0.6-28   stringi_1.4.3         RSQLite_2.1.1         pcaPP_1.9-73          foreach_1.4.4         checkmate_1.9.3       caTools_1.17.1.2      bibtex_0.4.2          Rdpack_0.11-0         SDMTools_1.1-221.1    rlang_0.3.4           pkgconfig_2.0.2      
##  [81] bitops_1.0-6          matrixStats_0.54.0    evaluate_0.14         lattice_0.20-38       ROCR_1.0-7            purrr_0.3.2           htmlwidgets_1.3       robust_0.4-18         cowplot_0.9.4         bit_1.1-14            tidyselect_0.2.5      plyr_1.8.4            magrittr_1.5          R6_2.4.0              fit.models_0.5-14     Hmisc_4.2-0           gplots_3.0.1.1        DBI_1.0.0             foreign_0.8-71        pillar_1.4.1         
## [101] withr_2.1.2           fitdistrplus_1.0-14   nnet_7.3-12           survival_2.44-1.1     tibble_2.1.3          future.apply_1.3.0    tsne_0.1-3            crayon_1.3.4          KernSmooth_2.23-15    plotly_4.9.0          rmarkdown_1.13        grid_3.6.0            data.table_1.12.2     blob_1.1.1            metap_1.1             digest_0.6.19         tidyr_0.8.3           R.utils_2.9.0         munsell_0.5.0         viridisLite_0.3.0
```
