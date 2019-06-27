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
##               17               20               24               24 
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
##   0           162           196          168
##   1            77           114          155
##   2           115           119           98
##   3            93            82           94
##   4            55            77           86
##   5            60            66           65
##   6            50            79           62
##   7            50            53           45
##   8            23            45           39
##   9            34            33           38
##   10           20            30           20
##   11           31            19           18
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
## Baiap2l1 7.352739e-235 1.1954540 0.714 0.014 9.419593e-231
## Cadps2   4.448365e-198 2.3426788 0.971 0.051 5.698801e-194
## Tbx3os2  1.753416e-182 0.6699683 0.443 0.005 2.246301e-178
## Cbln2    1.849806e-152 1.1941476 0.786 0.038 2.369787e-148
## Ntng1    1.302871e-150 1.0795688 0.686 0.028 1.669108e-146
## Ntrk2    1.157467e-148 2.0506117 0.971 0.074 1.482831e-144
```

```r
dim(markers)
```

```
## [1] 1474    5
```

```r
table(markers$avg_logFC > 0)
```

```
## 
## FALSE  TRUE 
##   688   786
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
## [1] 5941    7
```

```r
head(markers_all)
```

```
##                 p_val avg_logFC pct.1 pct.2     p_val_adj cluster    gene
## Tac1    2.713535e-256 1.8798988 0.885 0.384 3.476310e-252       0    Tac1
## Adcyap1 4.686385e-241 1.7593200 0.652 0.071 6.003728e-237       0 Adcyap1
## Celf4   2.211318e-236 1.5859987 0.876 0.376 2.832919e-232       0   Celf4
## Calca   2.095675e-228 1.3834990 0.934 0.365 2.684769e-224       0   Calca
## Syt4    1.738417e-195 1.2066331 0.944 0.771 2.227086e-191       0    Syt4
## Fgf13   2.568525e-195 0.9353625 0.973 0.834 3.290538e-191       0   Fgf13
```

```r
table(table(markers_all$gene))
```

```
## 
##    1    2    3    4    5    6    7 
## 1582  921  520  168   46    8    1
```

```r
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

dim(markers_all_single)
```

```
## [1] 1582    7
```

```r
table(table(markers_all_single$gene))
```

```
## 
##    1 
## 1582
```

```r
table(markers_all_single$cluster)
```

```
## 
##   0   1   3   4   5   6   7   8   9  10  11  12  13 
##  85 105  76 127 244 118 381  45  60 159  13  87  82
```

```r
head(markers_all_single)
```

```
##                 p_val avg_logFC pct.1 pct.2     p_val_adj cluster    gene
## Adcyap1 4.686385e-241  1.759320 0.652 0.071 6.003728e-237       0 Adcyap1
## Celf4   2.211318e-236  1.585999 0.876 0.376 2.832919e-232       0   Celf4
## Syt4    1.738417e-195  1.206633 0.944 0.771 2.227086e-191       0    Syt4
## Gal     6.083472e-195  1.779010 0.481 0.023 7.793537e-191       0     Gal
## Gpx3    2.671333e-152  1.437052 0.479 0.064 3.422245e-148       0    Gpx3
## Nrsn1   3.297842e-129  1.194130 0.783 0.450 4.224865e-125       0   Nrsn1
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
## slot for the RNA assay: Tusc2, mt-Nd5, Dusp15, Brap, Chd5, Syt4, Celf4
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
## Tac1    2.713535e-256 1.8798988 0.885 0.384 3.476310e-252       0    Tac1
## Adcyap1 4.686385e-241 1.7593200 0.652 0.071 6.003728e-237       0 Adcyap1
## Celf4   2.211318e-236 1.5859987 0.876 0.376 2.832919e-232       0   Celf4
## Calca   2.095675e-228 1.3834990 0.934 0.365 2.684769e-224       0   Calca
## Syt4    1.738417e-195 1.2066331 0.944 0.771 2.227086e-191       0    Syt4
## Fgf13   2.568525e-195 0.9353625 0.973 0.834 3.290538e-191       0   Fgf13
##         mean.in.cluster mean.out.of.cluster
## Tac1           2.749737           0.6976187
## Adcyap1        1.482089           0.1108635
## Celf4          2.362463           0.6189058
## Calca          3.214833           0.9258401
## Syt4           2.964189           1.6223466
## Fgf13          3.383159           2.1450409
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
## Car8          1.354856e-14  0.3956540 0.162 0.021 1.735705e-10
## Tmsb10        2.895653e-14  0.2941844 0.982 0.923 3.709621e-10
## Rpl21         1.471830e-10  0.3522771 0.921 0.744 1.885561e-06
## Rpl23a        2.644815e-10  0.2994179 0.931 0.759 3.388272e-06
## Rpl39         9.360848e-10  0.2654162 0.949 0.809 1.199218e-05
## Atp5e         5.415361e-09  0.3536987 0.675 0.458 6.937619e-05
## Pcp4          1.104481e-07  0.4766541 0.603 0.430 1.414951e-03
## S100a10       1.680102e-07  0.2816616 0.715 0.494 2.152379e-03
## Fxyd7         9.913820e-07  0.3208104 0.563 0.375 1.270059e-02
## Rps26         1.046363e-06  0.2727891 0.834 0.649 1.340496e-02
## Ndufa3        1.197237e-06  0.2601290 0.513 0.318 1.533781e-02
## Wdr89         1.749409e-06  0.2534736 0.238 0.108 2.241168e-02
## 1810058I24Rik 2.676280e-06  0.2817320 0.516 0.330 3.428582e-02
## Ndufv3        4.077302e-06  0.2586730 0.729 0.542 5.223432e-02
## Actb          6.570268e-06 -0.3934400 0.993 0.974 8.417170e-02
## Acbd5         1.371884e-05  0.2559061 0.390 0.244 1.757521e-01
## S100a11       2.872269e-05  0.2907921 0.477 0.320 3.679664e-01
## Wdfy1         5.135154e-04 -0.3368053 0.072 0.155 1.000000e+00
## Polr2f        5.611144e-04  0.2986471 0.235 0.143 1.000000e+00
## Llph          5.811095e-04  0.2524122 0.184 0.100 1.000000e+00
## Gpx3          5.920469e-04  0.2549094 0.574 0.434 1.000000e+00
## Kansl1        1.013676e-03  0.2540896 0.166 0.090 1.000000e+00
## Necab1        1.068096e-03 -0.4008516 0.195 0.284 1.000000e+00
## Chchd10       1.804446e-03  0.2886958 0.170 0.096 1.000000e+00
## Arl6ip1       2.382005e-03 -0.2513330 0.700 0.688 1.000000e+00
## Arf3          1.021658e-02 -0.3268018 0.549 0.552 1.000000e+00
## Ptgir         2.212238e-02 -0.2652934 0.072 0.120 1.000000e+00
## Rnf146        3.397373e-02 -0.2669655 0.079 0.126 1.000000e+00
## Etv5          4.018323e-02 -0.3213612 0.173 0.222 1.000000e+00
## Vdac1         5.554814e-02 -0.2996588 0.415 0.425 1.000000e+00
## Ist1          6.127142e-02 -0.2639362 0.094 0.136 1.000000e+00
## Nudcd3        6.745731e-02 -0.2508401 0.112 0.155 1.000000e+00
## Aplp2         6.773344e-02 -0.2858165 0.606 0.563 1.000000e+00
## Cbx3          7.385992e-02 -0.2534048 0.347 0.377 1.000000e+00
## Aff4          9.446635e-02 -0.2550704 0.195 0.236 1.000000e+00
## Usp22         9.659045e-02 -0.2980711 0.289 0.313 1.000000e+00
## Nacc2         1.316943e-01 -0.2985409 0.350 0.356 1.000000e+00
## Amer2         1.371509e-01 -0.2645151 0.159 0.189 1.000000e+00
## Gpsm3         1.580521e-01 -0.2550681 0.119 0.148 1.000000e+00
## Eif5          1.593214e-01 -0.2566491 0.408 0.406 1.000000e+00
## Rapgef4       1.723811e-01 -0.2626167 0.170 0.198 1.000000e+00
## Bhlhe41       1.891923e-01 -0.2806677 0.415 0.423 1.000000e+00
## Gng5          1.988570e-01 -0.2598815 0.177 0.205 1.000000e+00
## Alkbh5        2.042601e-01 -0.3241269 0.181 0.207 1.000000e+00
## Hdlbp         2.138709e-01 -0.2762089 0.097 0.120 1.000000e+00
## Tsc22d3       2.141794e-01 -0.2733712 0.271 0.293 1.000000e+00
## Clasp2        2.637078e-01 -0.2660012 0.188 0.207 1.000000e+00
## Sult4a1       3.116372e-01 -0.2618857 0.531 0.487 1.000000e+00
## Nrn1          3.475281e-01 -0.2559416 0.412 0.396 1.000000e+00
## Tspan2        4.079501e-01 -0.2546101 0.184 0.193 1.000000e+00
## Abhd2         5.065008e-01 -0.3497591 0.412 0.392 1.000000e+00
## Setd3         5.437484e-01 -0.2677854 0.242 0.238 1.000000e+00
## Spry2         7.045794e-01 -0.2938026 0.191 0.191 1.000000e+00
## Pam           7.095311e-01 -0.5316723 0.581 0.523 1.000000e+00
## Rgs4          8.191299e-01 -0.2721859 0.458 0.427 1.000000e+00
## Nuak1         8.281399e-01 -0.2804775 0.173 0.167 1.000000e+00
## Atp1b1        9.304997e-01 -0.2569386 0.357 0.317 1.000000e+00
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
## assay: Nuak1, Setd3, Tspan2, Sult4a1, Clasp2, Tsc22d3, Hdlbp, Gng5, Eif5,
## Amer2, Usp22, Aff4, Nudcd3, Ist1, Vdac1, Rnf146, Arf3, Arl6ip1, Necab1,
## Kansl1, Polr2f, Wdfy1, Acbd5, 1810058I24Rik, Wdr89, Ndufa3, Rps26, Atp5e,
## Rpl39, Rpl23a, Rpl21, Tmsb10
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
