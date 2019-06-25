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
## Xkr4              0.003732825             0.004357195            0.007326931
## Sox17             0.000000000             0.014117159            0.000000000
## Mrpl15            0.161265269             0.125769122            0.172108012
## Lypla1            0.185685662             0.156808334            0.240187355
## Tcea1             0.197782618             0.234938921            0.229517542
## Rgs20             0.062002520             0.100209702            0.078421134
```

```r
contr <- makeContrasts(orig.identUCD_Supp_VitE - orig.identUCD_Adj_VitE, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrasts = contr)
tmp <- eBayes(tmp)
topTable(tmp, sort.by = "P", n = 20) # top 20 DE genes
```

```
##               logFC   AveExpr         t      P.Value    adj.P.Val         B
## Rpl21    -0.6320258 1.9936203 -6.500350 1.420406e-10 1.789995e-06 13.614158
## Rpl17    -0.5422245 1.8828905 -5.719612 1.515483e-08 6.529146e-05  9.207235
## Rpl23a   -0.5156828 2.1544077 -5.715133 1.554312e-08 6.529146e-05  9.183435
## Pcp4     -0.6955342 1.3863141 -5.560079 3.692219e-08 1.163234e-04  8.370104
## Rpl24    -0.4706538 2.0936917 -5.344913 1.185232e-07 2.987258e-04  7.275693
## H3f3b    -0.5030337 1.9814042 -5.208081 2.436998e-07 4.917478e-04  6.600563
## Rpl39    -0.4432437 2.1927211 -5.186137 2.731499e-07 4.917478e-04  6.493809
## Rps8     -0.4599371 2.4562099 -5.094602 4.376310e-07 6.893782e-04  6.053035
## Tmsb10   -0.3874920 3.2659606 -5.003088 6.959682e-07 9.567369e-04  5.619698
## Rpl32    -0.4505451 2.3814281 -4.985773 7.591945e-07 9.567369e-04  5.538535
## Rpl23    -0.4276946 2.0730354 -4.897481 1.177984e-06 1.349542e-03  5.128776
## Hspa8    -0.4598663 2.2022393 -4.704895 2.998937e-06 3.149383e-03  4.258889
## Fau      -0.3608879 2.7663506 -4.673175 3.486960e-06 3.380206e-03  4.118771
## mt-Co3   -0.3499016 4.1954016 -4.646346 3.958495e-06 3.415486e-03  4.000955
## Map1lc3a -0.3896356 2.7386762 -4.632454 4.226139e-06 3.415486e-03  3.940201
## Rpl36    -0.4136916 2.1184903 -4.620437 4.471586e-06 3.415486e-03  3.887789
## Rps24    -0.4223393 1.5455986 -4.614055 4.607464e-06 3.415486e-03  3.860002
## Rps15    -0.4227129 1.9978401 -4.597629 4.975613e-06 3.483482e-03  3.788658
## Rpl6     -0.4119715 2.2377053 -4.585744 5.259397e-06 3.488365e-03  3.737191
## Atp6v1c1 -0.3549174 0.7340975 -4.527828 6.879465e-06 4.176513e-03  3.488179
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
## 		 the algorithm is scoring 3126 nontrivial nodes
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
## 	 Level 15:	16 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 14:	31 nodes to be scored	(36 eliminated genes)
```

```
## 
## 	 Level 13:	84 nodes to be scored	(43 eliminated genes)
```

```
## 
## 	 Level 12:	130 nodes to be scored	(427 eliminated genes)
```

```
## 
## 	 Level 11:	203 nodes to be scored	(541 eliminated genes)
```

```
## 
## 	 Level 10:	289 nodes to be scored	(582 eliminated genes)
```

```
## 
## 	 Level 9:	357 nodes to be scored	(1209 eliminated genes)
```

```
## 
## 	 Level 8:	407 nodes to be scored	(1604 eliminated genes)
```

```
## 
## 	 Level 7:	487 nodes to be scored	(1819 eliminated genes)
```

```
## 
## 	 Level 6:	475 nodes to be scored	(2298 eliminated genes)
```

```
## 
## 	 Level 5:	334 nodes to be scored	(2399 eliminated genes)
```

```
## 
## 	 Level 4:	195 nodes to be scored	(2479 eliminated genes)
```

```
## 
## 	 Level 3:	88 nodes to be scored	(2827 eliminated genes)
```

```
## 
## 	 Level 2:	19 nodes to be scored	(2827 eliminated genes)
```

```
## 
## 	 Level 1:	1 nodes to be scored	(2827 eliminated genes)
```

```r
	GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
```

```
##         GO.ID                                                            Term Annotated Significant Expected  Fisher
## 1  GO:0006412                                                     translation       517          56     7.32 3.0e-23
## 2  GO:0002181                                         cytoplasmic translation        77          21     1.09 8.2e-22
## 3  GO:0000028                                ribosomal small subunit assembly        18           9     0.25 8.0e-13
## 4  GO:0000027                                ribosomal large subunit assembly        29           8     0.41 4.5e-09
## 5  GO:1902255 positive regulation of intrinsic apoptotic signaling pathway...         6           4     0.08 5.7e-07
## 6  GO:1904667        negative regulation of ubiquitin protein ligase activity         7           3     0.10 9.4e-05
## 7  GO:0047497                       mitochondrion transport along microtubule        19           4     0.27 0.00013
## 8  GO:0071353                              cellular response to interleukin-4        21           4     0.30 0.00019
## 9  GO:0097214          positive regulation of lysosomal membrane permeability         2           2     0.03 0.00020
## 10 GO:0007409                                                    axonogenesis       353          14     5.00 0.00049
## 11 GO:0006880                          intracellular sequestering of iron ion         3           2     0.04 0.00059
## 12 GO:2000582 positive regulation of ATP-dependent microtubule motor activ...         3           2     0.04 0.00059
## 13 GO:0000462 maturation of SSU-rRNA from tricistronic rRNA transcript (SS...        30           4     0.42 0.00080
## 14 GO:0006364                                                 rRNA processing       177          12     2.51 0.00100
## 15 GO:0050848                        regulation of calcium-mediated signaling        60           5     0.85 0.00156
## 16 GO:0061844 antimicrobial humoral immune response mediated by antimicrob...        17           3     0.24 0.00164
## 17 GO:0006874                                cellular calcium ion homeostasis       272          11     3.85 0.00168
## 18 GO:0002227                                innate immune response in mucosa         5           2     0.07 0.00194
## 19 GO:0071635 negative regulation of transforming growth factor beta produ...         5           2     0.07 0.00194
## 20 GO:0032272                   negative regulation of protein polymerization        63           5     0.89 0.00194
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



