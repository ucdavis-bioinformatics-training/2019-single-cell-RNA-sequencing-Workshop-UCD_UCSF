---
title: "Single Cell RNAseq Part 5"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

## Load libraries

```r
library(Seurat)
library(ggplot2)
```

## Load the Seurat object

```r
load(file="pca_sample_corrected.RData")
experiment.aggregate
```

```
## An object of class Seurat 
## 12811 features across 2681 samples within 1 assay 
## Active assay: RNA (12811 features)
##  1 dimensional reduction calculated: pca
```

## Identifying clusters

Seurat implements an graph-based clustering approach. Distances between the cells are calculated based on previously identified PCs. Seurat approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNAseq data. Briefly, Seurat identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors (KNN) and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B.

The FindClusters function implements the procedure, and contains a resolution parameter that sets the granularity of the downstream clustering, with increased values leading to a greater number of clusters. I tend to like to perform a series of resolutions, investigate and choose.


```r
use.pcs = 1:29 

?FindNeighbors
experiment.aggregate <- FindNeighbors(experiment.aggregate, reduction="pca", dims = use.pcs)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
?FindCluster
```

```
## No documentation for 'FindCluster' in specified packages and libraries:
## you could try '??FindCluster'
```

```r
experiment.aggregate <- FindClusters(
    object = experiment.aggregate, 
    resolution = seq(0.25,4,0.25), 
    verbose = FALSE
)
```

Lets first investigate how many clusters each resolution produces and set it to the smallest resolutions of 0.5 (fewest clusters). 


```r
sapply(grep("res",colnames(experiment.aggregate@meta.data),value = TRUE),
       function(x) length(unique(experiment.aggregate@meta.data[,x])))
```

```
## RNA_snn_res.0.25  RNA_snn_res.0.5 RNA_snn_res.0.75    RNA_snn_res.1 
##               10               14               15               16 
## RNA_snn_res.1.25  RNA_snn_res.1.5 RNA_snn_res.1.75    RNA_snn_res.2 
##               17               21               24               24 
## RNA_snn_res.2.25  RNA_snn_res.2.5 RNA_snn_res.2.75    RNA_snn_res.3 
##               25               26               27               28 
## RNA_snn_res.3.25  RNA_snn_res.3.5 RNA_snn_res.3.75    RNA_snn_res.4 
##               28               28               28               29
```

```r
Idents(experiment.aggregate) <- "RNA_snn_res.0.5"
```

Finally,  lets produce a table of cluster to sample assignments.

```r
table(Idents(experiment.aggregate),experiment.aggregate$orig.ident)
```

```
##     
##      UCD_Adj_VitE UCD_Supp_VitE UCD_VitE_Def
##   0           155           173          158
##   1            78           116          156
##   2           115           123          105
##   3            93            81           92
##   4            57           101           67
##   5            55            77           86
##   6            60            66           65
##   7            48            49           45
##   8            23            45           39
##   9            35            33           37
##   10           20            30           21
##   11           31            19           17
##   12           18            20           19
##   13           20            14           19
```

tSNE dimensionality reduction plots are then used to visualise clustering results. As input to the tSNE, you should use the same PCs as input to the clustering analysis.


```r
experiment.aggregate <- RunTSNE(
  object = experiment.aggregate,
  reduction.use = "pca",
  dims.use = use.pcs,
  do.fast = TRUE)
```

Plot TSNE coloring by the slot 'ident' (default).

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


Plot TSNE coloring by the clustering resolution 4

```r
DimPlot(object = experiment.aggregate, group.by="RNA_snn_res.4", pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

FeaturePlot can be used to color cells with a 'feature', non categorical data, like number of UMIs

```r
FeaturePlot(experiment.aggregate, features = c('nCount_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-9-1.png)<!-- -->
and number of genes present

```r
FeaturePlot(experiment.aggregate, features = c('nFeature_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

percent mitochondrial 

```r
FeaturePlot(experiment.aggregate, features = c('percent.mito'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

TSNE plot by cell cycle

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, group.by = "cell.cycle", reduction = "tsne" )
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


## Building  a  tree relating the 'average' cell from each cluster. Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space.


```r
experiment.aggregate <- BuildClusterTree(
  experiment.aggregate, dims = use.pcs)

PlotClusterTree(experiment.aggregate)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


```r
DimPlot(object = experiment.aggregate, pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


```r
experiment.merged <- RenameIdents(
  object = experiment.aggregate,
  '2' = '0'
)
DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

Plot TSNE coloring by the slot 'orig.ident' (sample names) with alpha colors turned on.

```r
DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "tsne" )
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```r
## Pretty tsne using alpha
p <- DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "tsne", do.return = T)
alpha.use <- 2/5
p$layers[[1]]$mapping$alpha <- alpha.use
p + scale_alpha_continuous(range = alpha.use, guide = F)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-16-2.png)<!-- -->

## Identifying Marker Genes

Seurat can help you find markers that define clusters via differential expression.

`FindMarkers` identifies markers for a cluster relative to all other clusters.

`FindAllMarkers` does so for all clusters

`FindAllMarkersNode` defines all markers that split a Node __(Warning: need to validate)__


```r
?FindMarkers

markers = FindMarkers(experiment.merged, ident.1=c(10), genes.use = VariableFeatures(experiment.merged))

head(markers)
```

```
##                  p_val avg_logFC pct.1 pct.2     p_val_adj
## Baiap2l1 2.366963e-241 1.2648290 0.718 0.014 3.032316e-237
## Cadps2   2.382879e-201 2.3820185 0.972 0.051 3.052706e-197
## Tbx3os2  8.086208e-180 0.6629623 0.437 0.005 1.035924e-175
## Cbln2    4.633891e-150 1.1834309 0.775 0.038 5.936478e-146
## Ntng1    2.787527e-148 1.0695523 0.676 0.028 3.571101e-144
## Ntrk2    3.276192e-146 2.0375267 0.958 0.074 4.197130e-142
```

```r
dim(markers)
```

```
## [1] 1499    5
```

```r
table(markers$avg_logFC > 0)
```

```
## 
## FALSE  TRUE 
##   698   801
```

 
pct.1 and pct.2 are the proportion of cells with expression above 0 in ident.1 and ident.2 respectively. p_val is the raw p_value associated with the differntial expression test with adjusted value in p_val_adj. avg_logFC is the average log fold change difference between the two groups. 
 
avg_diff (lines 130, 193 and) appears to be the difference in log(x = mean(x = exp(x = x) - 1) + 1) between groups.  It doesn’t seem like this should work out to be the signed ratio of pct.1 to pct.2 so I must be missing something.  It doesn’t seem to be related at all to how the p-values are calculated so maybe it doesn’t matter so much, and the sign is probably going to be pretty robust to how expression is measured.

Can use a violin plot to visualize the expression pattern of some markers

```r
VlnPlot(object = experiment.merged, features = rownames(markers)[1:2], pt.size = 0.05)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

Or a feature plot

```r
FeaturePlot(
    experiment.merged, 
    head(rownames(markers), n=6), 
    cols = c("lightgrey", "blue"), 
    ncol = 2
)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

```r
FeaturePlot(    
    experiment.merged, 
    "Fxyd1", 
    cols = c("lightgrey", "blue") 
)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

FindAllMarkers can be used to automate the process across all genes.
__WARNING: TAKES A LONG TIME TO RUN__


```r
markers_all <- FindAllMarkers(
    object = experiment.merged, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    thresh.use = 0.25
)
```

```
## Calculating cluster 0
```

```
## Calculating cluster 1
```

```
## Calculating cluster 3
```

```
## Calculating cluster 4
```

```
## Calculating cluster 5
```

```
## Calculating cluster 6
```

```
## Calculating cluster 7
```

```
## Calculating cluster 8
```

```
## Calculating cluster 9
```

```
## Calculating cluster 10
```

```
## Calculating cluster 11
```

```
## Calculating cluster 12
```

```
## Calculating cluster 13
```

```r
dim(markers_all)
```

```
## [1] 6027    7
```

```r
head(markers_all)
```

```
##                 p_val avg_logFC pct.1 pct.2     p_val_adj cluster    gene
## Tac1    5.714341e-268 1.8979687 0.905 0.383 7.320642e-264       0    Tac1
## Adcyap1 2.016326e-252 1.7727603 0.672 0.071 2.583116e-248       0 Adcyap1
## Celf4   5.793510e-243 1.5786648 0.894 0.376 7.422065e-239       0   Celf4
## Calca   3.802333e-222 1.3510992 0.941 0.371 4.871169e-218       0   Calca
## Gal     2.792374e-204 1.8042863 0.497 0.023 3.577310e-200       0     Gal
## Fgf13   7.264395e-199 0.9410338 0.982 0.833 9.306417e-195       0   Fgf13
```

```r
table(table(markers_all$gene))
```

```
## 
##    1    2    3    4    5    6    7 
## 1596  932  538  162   50    8    1
```

```r
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

dim(markers_all_single)
```

```
## [1] 1596    7
```

```r
table(table(markers_all_single$gene))
```

```
## 
##    1 
## 1596
```

```r
table(markers_all_single$cluster)
```

```
## 
##   0   1   3   4   5   6   7   8   9  10  11  12  13 
##  91 103  81 133 126 244 365  44  59 163  14  89  84
```

```r
head(markers_all_single)
```

```
##                  p_val avg_logFC pct.1 pct.2     p_val_adj cluster
## Adcyap1  2.016326e-252  1.772760 0.672 0.071 2.583116e-248       0
## Gal      2.792374e-204  1.804286 0.497 0.023 3.577310e-200       0
## Syt4     1.693291e-197  1.198862 0.954 0.769 2.169275e-193       0
## Gpx3     7.901586e-161  1.464333 0.495 0.063 1.012272e-156       0
## Nrsn1    7.388698e-141  1.225448 0.806 0.445 9.465661e-137       0
## Cacna2d1 2.530465e-104  1.049548 0.749 0.491 3.241779e-100       0
##              gene
## Adcyap1   Adcyap1
## Gal           Gal
## Syt4         Syt4
## Gpx3         Gpx3
## Nrsn1       Nrsn1
## Cacna2d1 Cacna2d1
```

Plot a heatmap of genes by cluster for the top 5 marker genes per cluster

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
top5 <- markers_all_single %>% group_by(cluster) %>% top_n(5, avg_logFC)
dim(top5)
```

```
## [1] 65  7
```

```r
DoHeatmap(
    object = experiment.merged, 
    features = top5$gene
) 
```

```
## Warning in DoHeatmap(object = experiment.merged, features = top5$gene): The
## following features were omitted as they were not found in the scale.data
## slot for the RNA assay: Tusc2, Dusp15, Brap, Chd5, Nrsn1
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-21-1.png)<!-- -->


```r
# Get expression of genes for cells in and out of each cluster
getGeneClusterMeans <- function(gene, cluster){
  x <- GetAssayData(experiment.merged)[gene,]
  m <- tapply(x, ifelse(Idents(experiment.merged) == cluster, 1, 0), mean)
  mean.in.cluster <- m[2]
  mean.out.of.cluster <- m[1]
  return(list(mean.in.cluster = mean.in.cluster, mean.out.of.cluster = mean.out.of.cluster))
}

## for sake of time only using first six (head)
means <- mapply(getGeneClusterMeans, head(markers_all[,"gene"]), head(markers_all[,"cluster"]))
means <- matrix(unlist(means), ncol = 2, byrow = T)

colnames(means) <- c("mean.in.cluster", "mean.out.of.cluster")
rownames(means) <- head(markers_all[,"gene"])
markers_all2 <- cbind(head(markers_all), means)
head(markers_all2)
```

```
##                 p_val avg_logFC pct.1 pct.2     p_val_adj cluster    gene
## Tac1    5.714341e-268 1.8979687 0.905 0.383 7.320642e-264       0    Tac1
## Adcyap1 2.016326e-252 1.7727603 0.672 0.071 2.583116e-248       0 Adcyap1
## Celf4   5.793510e-243 1.5786648 0.894 0.376 7.422065e-239       0   Celf4
## Calca   3.802333e-222 1.3510992 0.941 0.371 4.871169e-218       0   Calca
## Gal     2.792374e-204 1.8042863 0.497 0.023 3.577310e-200       0     Gal
## Fgf13   7.264395e-199 0.9410338 0.982 0.833 9.306417e-195       0   Fgf13
##         mean.in.cluster mean.out.of.cluster
## Tac1           2.816402          0.69991158
## Adcyap1        1.526563          0.11242747
## Celf4          2.407467          0.62606306
## Calca          3.233145          0.95348629
## Gal            1.150708          0.03477168
## Fgf13          3.419593          2.14811956
```

## Finishing up clusters.

At this point in time you should use the tree, markers, domain knowledge, and goals to finalize your clusters. This may mean adjusting PCA to use, mergers clusters together, choosing a new resolutions, etc. When finished you can further name it cluster by something more informative. Ex.

```r
experiment.clusters <- experiment.merged
experiment.clusters <- RenameIdents(
  object = experiment.clusters,
  '0' = 'cell_type_A'
)

DimPlot(object = experiment.clusters, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
experiment.merged$finalcluster <- Idents(experiment.merged)
```

## Subsetting samples

```r
experiment.sample2 <- subset(experiment.merged, orig.ident == "UCD_Supp_VitE")

DimPlot(object = experiment.sample2, group.by = "RNA_snn_res.0.5", pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

```r
FeaturePlot(experiment.sample2, features =c('Calca'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-24-2.png)<!-- -->

```r
FeaturePlot(experiment.sample2, features =c('Adcyap1'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-24-3.png)<!-- -->

### Adding in a new metadata column representing samples within clusters


```r
experiment.merged$samplecluster = paste(experiment.merged$orig.ident,experiment.merged$finalcluster,sep = '-')

# set the identity to the new variable 
Idents(experiment.merged) <- "samplecluster"

markers.comp <- FindMarkers(experiment.merged, ident.1 = "UCD_Adj_VitE-0", ident.2= c("UCD_Supp_VitE-0","UCD_VitE_Def-0"))

markers.comp
```

```
##                      p_val  avg_logFC pct.1 pct.2    p_val_adj
## Car8          6.108801e-14  0.3874529 0.163 0.021 7.825984e-10
## Tmsb10        2.584463e-13  0.2793655 0.985 0.934 3.310956e-09
## Rpl39         9.980978e-10  0.2625274 0.963 0.823 1.278663e-05
## Rpl23a        1.635937e-09  0.2771557 0.948 0.782 2.095798e-05
## Rpl21         2.019520e-09  0.3275539 0.937 0.775 2.587206e-05
## Atp5e         1.350918e-08  0.3454480 0.685 0.470 1.730660e-04
## Pcp4          1.742560e-07  0.4670440 0.615 0.442 2.232394e-03
## Fxyd7         2.521184e-07  0.3478833 0.574 0.376 3.229889e-03
## S100a10       3.765650e-07  0.2693872 0.730 0.510 4.824175e-03
## Ahsa2         1.149068e-06  0.2569993 0.274 0.129 1.472071e-02
## Ndufa3        1.532781e-06  0.2508090 0.522 0.326 1.963646e-02
## Wdr89         2.577710e-06  0.2522381 0.244 0.113 3.302304e-02
## Rps26         3.809858e-06  0.2594906 0.848 0.673 4.880809e-02
## 1810058I24Rik 5.495249e-06  0.2733952 0.522 0.338 7.039964e-02
## Acbd5         1.467241e-05  0.2575561 0.389 0.243 1.879683e-01
## Actb          3.905710e-05 -0.3837563 0.996 0.975 5.003605e-01
## S100a11       4.898904e-05  0.2824516 0.489 0.333 6.275986e-01
## Necab1        5.999469e-04 -0.4273210 0.196 0.292 1.000000e+00
## Wdfy1         7.853582e-04 -0.3043508 0.074 0.156 1.000000e+00
## Kansl1        1.056426e-03  0.2649178 0.167 0.089 1.000000e+00
## Polr2f        1.208560e-03  0.2707323 0.237 0.148 1.000000e+00
## Chchd10       1.619698e-03  0.3004953 0.174 0.098 1.000000e+00
## Erp29         3.306382e-03  0.2684310 0.274 0.186 1.000000e+00
## Arf3          4.256416e-03 -0.3421780 0.537 0.553 1.000000e+00
## Ptgir         1.345200e-02 -0.2951010 0.070 0.123 1.000000e+00
## Plp1          1.792279e-02 -0.4266177 0.119 0.175 1.000000e+00
## Cadm3         4.880916e-02 -0.2508777 0.133 0.181 1.000000e+00
## Etv5          4.909260e-02 -0.3095786 0.178 0.225 1.000000e+00
## Cbx3          5.095870e-02 -0.2695332 0.352 0.388 1.000000e+00
## Nudcd3        6.616364e-02 -0.2527923 0.115 0.159 1.000000e+00
## Aplp2         6.826221e-02 -0.2766304 0.607 0.562 1.000000e+00
## Vdac1         8.265813e-02 -0.2782682 0.419 0.424 1.000000e+00
## Aff4          9.796625e-02 -0.2534300 0.200 0.242 1.000000e+00
## Vim           1.008757e-01 -0.3129994 0.156 0.193 1.000000e+00
## Gpsm3         1.011990e-01 -0.2799424 0.119 0.154 1.000000e+00
## Gng5          1.426415e-01 -0.2781454 0.178 0.211 1.000000e+00
## Amer2         1.471798e-01 -0.2571606 0.163 0.193 1.000000e+00
## Alkbh5        1.920259e-01 -0.3234712 0.181 0.209 1.000000e+00
## Clasp2        2.003888e-01 -0.2811156 0.193 0.216 1.000000e+00
## Nacc2         2.181322e-01 -0.2625106 0.356 0.354 1.000000e+00
## Evi5          2.409357e-01 -0.2514888 0.189 0.209 1.000000e+00
## Bhlhe41       2.659677e-01 -0.2678791 0.426 0.428 1.000000e+00
## Pam           5.178273e-01 -0.5162452 0.589 0.540 1.000000e+00
## Etv1          6.145764e-01 -0.2639596 0.163 0.168 1.000000e+00
## Spry2         6.310345e-01 -0.2995632 0.193 0.195 1.000000e+00
## Abhd2         6.328808e-01 -0.3177709 0.422 0.397 1.000000e+00
## Rgs4          8.109231e-01 -0.3109885 0.452 0.422 1.000000e+00
```

```r
DoHeatmap(experiment.merged,
          cells = rownames(experiment.merged@meta.data)[experiment.merged@meta.data$samplecluster %in% c( "UCD_Adj_VitE-0", "UCD_Supp_VitE-0" )],
          features = rownames(markers.comp),
          )
```

```
## Warning in DoHeatmap(experiment.merged,
## cells = rownames(experiment.merged@meta.data)
## [experiment.merged@meta.data$samplecluster %in% : The following features
## were omitted as they were not found in the scale.data slot for the RNA
## assay: Evi5, Clasp2, Amer2, Gng5, Vim, Aff4, Vdac1, Nudcd3, Arf3, Erp29,
## Polr2f, Kansl1, Wdfy1, Necab1, Acbd5, 1810058I24Rik, Rps26, Wdr89, Ndufa3,
## Ahsa2, Atp5e, Rpl21, Rpl23a, Rpl39, Tmsb10
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

```r
Idents(experiment.merged) <- "finalcluster"
```

And last lets save all the objects in our session.

```r
save(list=ls(), file="clusters_seurat_object.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/master/scrnaseq_analysis/scRNA_Workshop-PART6.Rmd", "scRNA_Workshop-PART6.Rmd")
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
## [1] dplyr_0.8.1   ggplot2_3.2.0 Seurat_3.0.2 
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
## [37] withr_2.1.2         ROCR_1.0-7          pbapply_1.4-0      
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
## [79] future.apply_1.3.0  KernSmooth_2.23-15  metap_1.1          
## [82] ape_5.3             stringi_1.4.3       parallel_3.6.0     
## [85] Rcpp_1.0.1          sctransform_0.2.0   png_0.1-7          
## [88] tidyselect_0.2.5    xfun_0.7            lmtest_0.9-37
```
