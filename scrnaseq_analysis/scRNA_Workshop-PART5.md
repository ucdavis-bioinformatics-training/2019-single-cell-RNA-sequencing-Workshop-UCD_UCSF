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
##               18               21               24               24 
## RNA_snn_res.2.25  RNA_snn_res.2.5 RNA_snn_res.2.75    RNA_snn_res.3 
##               25               26               26               27 
## RNA_snn_res.3.25  RNA_snn_res.3.5 RNA_snn_res.3.75    RNA_snn_res.4 
##               28               27               28               29
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
##   0           157           174          161
##   1           118           141          106
##   2            77           115          155
##   3            94            82           94
##   4            55            77           86
##   5            51            79           62
##   6            60            66           65
##   7            50            52           45
##   8            23            45           39
##   9            34            33           37
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

Merge Clustering results

```r
experiment.merged = experiment.aggregate
# originally set clusters to resolutionm 0.5
Idents(experiment.merged) <- "RNA_snn_res.0.5"

table(Idents(experiment.merged))
```

```
## 
##   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
## 492 365 347 270 218 192 191 147 107 104  70  68  57  53
```

```r
# based on TSNE and Heirarchical tree
# merge clusters 6 and 7 into 0 and cluster 9 into 13
experiment.merged <- RenameIdents(
  object = experiment.merged,
  '6' = '0', '7' = '0', '9' = '13'
)

table(Idents(experiment.merged))
```

```
## 
##   0  13   1   2   3   4   5   8  10  11  12 
## 830 157 365 347 270 218 192 107  70  68  57
```

```r
DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
experiment.examples <- experiment.merged
# in order to reporder the clusters for plotting purposes
# take a look at the levels, which indicates the ordering
levels(experiment.examples@active.ident)
```

```
##  [1] "0"  "13" "1"  "2"  "3"  "4"  "5"  "8"  "10" "11" "12"
```

```r
# relevel setting 5 to the first factor
experiment.examples@active.ident <- relevel(experiment.examples@active.ident, "5")
levels(experiment.examples@active.ident)
```

```
##  [1] "5"  "0"  "13" "1"  "2"  "3"  "4"  "8"  "10" "11" "12"
```

```r
# now cluster 5 is the "first" factor

# relevel all the factors to the order I want
Idents(experiment.examples) <- factor(experiment.examples@active.ident, levels=c("5","13","1","2","3","0","4","8","11","12","10","14"))
levels(experiment.examples@active.ident)
```

```
##  [1] "5"  "13" "1"  "2"  "3"  "0"  "4"  "8"  "11" "12" "10"
```

```r
DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-15-2.png)<!-- -->

```r
### Re-assign clustering result to resolution 4 for cells in cluster 0 (@ reslution 0.5) [adding a R prefix]
newIdent = as.character(Idents(experiment.examples))
newIdent[newIdent == '0'] = paste0("R",as.character(experiment.examples$RNA_snn_res.4[newIdent == '0']))

Idents(experiment.examples) <- as.factor(newIdent)

table(Idents(experiment.examples))
```

```
## 
##   1  10  11  12  13   2   3   4   5   8  R1 R10 R12 R13 R15 R16 R21 R26 
## 365  70  68  57 157 347 270 218 192 107 180   1  87  83  73   3  58   1 
## R27 R28  R5  R7  R8  R9 
##  35   1 151 118  37   2
```

```r
DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-15-3.png)<!-- -->

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

Removing cells assigned to clusters from a plot, So here plot all clusters but clusters "3" and "5"

```r
# create a new tmp object with those removed 
experiment.aggregate.tmp <- experiment.aggregate[,-which(Idents(experiment.aggregate) %in% c("3","5"))]

dim(experiment.aggregate)
```

```
## [1] 12811  2681
```

```r
dim(experiment.aggregate.tmp)
```

```
## [1] 12811  2219
```

```r
DimPlot(object = experiment.aggregate.tmp, group.by="orig.ident", pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

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

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

Or a feature plot

```r
FeaturePlot(
    experiment.merged, 
    head(rownames(markers), n=6), 
    cols = c("lightgrey", "blue"), 
    ncol = 2
)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
FeaturePlot(    
    experiment.merged, 
    "Fxyd1", 
    cols = c("lightgrey", "blue") 
)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-20-2.png)<!-- -->

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
## Calculating cluster 13
```

```
## Calculating cluster 1
```

```
## Calculating cluster 2
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
## Calculating cluster 8
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

```r
dim(markers_all)
```

```
## [1] 4140    7
```

```r
head(markers_all)
```

```
##                 p_val avg_logFC pct.1 pct.2     p_val_adj cluster    gene
## Pcp4    1.850077e-112 0.9957707 0.630 0.197 2.370134e-108       0    Pcp4
## Tac1     3.148614e-60 0.6683619 0.737 0.458  4.033689e-56       0    Tac1
## Marcks   3.243647e-59 0.7666353 0.681 0.397  4.155436e-55       0  Marcks
## Adcyap1  3.955847e-53 0.9475240 0.433 0.178  5.067836e-49       0 Adcyap1
## Nrsn1    2.536293e-50 0.8374397 0.716 0.486  3.249245e-46       0   Nrsn1
## Gal      1.310030e-48 0.8974830 0.324 0.100  1.678279e-44       0     Gal
```

```r
table(table(markers_all$gene))
```

```
## 
##    1    2    3    4    5    6 
## 1511  689  274   87   15    1
```

```r
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

dim(markers_all_single)
```

```
## [1] 1511    7
```

```r
table(table(markers_all_single$gene))
```

```
## 
##    1 
## 1511
```

```r
table(markers_all_single$cluster)
```

```
## 
##   0  13   1   2   3   4   5   8  10  11  12 
##  19 112  98 232 148 150 224  55 338  20 115
```

```r
head(markers_all_single)
```

```
##                p_val avg_logFC pct.1 pct.2    p_val_adj cluster    gene
## Nrsn1   2.536293e-50 0.8374397 0.716 0.486 3.249245e-46       0   Nrsn1
## Gm13889 1.365278e-37 0.6624093 0.658 0.491 1.749057e-33       0 Gm13889
## Ctnnd2  4.716368e-20 0.5710907 0.339 0.201 6.042138e-16       0  Ctnnd2
## Cpeb2   3.676886e-11 0.4783811 0.500 0.444 4.710458e-07       0   Cpeb2
## Nrip1   4.853869e-10 0.6027342 0.290 0.217 6.218291e-06       0   Nrip1
## Gprasp1 7.501332e-10 0.4571512 0.408 0.328 9.609956e-06       0 Gprasp1
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
## [1] 55  7
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
## slot for the RNA assay: Mest, Nwd2, Pik3r1, Nrip1, Cpeb2, Ctnnd2, Gm13889,
## Nrsn1
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-22-1.png)<!-- -->


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
## Pcp4    1.850077e-112 0.9957707 0.630 0.197 2.370134e-108       0    Pcp4
## Tac1     3.148614e-60 0.6683619 0.737 0.458  4.033689e-56       0    Tac1
## Marcks   3.243647e-59 0.7666353 0.681 0.397  4.155436e-55       0  Marcks
## Adcyap1  3.955847e-53 0.9475240 0.433 0.178  5.067836e-49       0 Adcyap1
## Nrsn1    2.536293e-50 0.8374397 0.716 0.486  3.249245e-46       0   Nrsn1
## Gal      1.310030e-48 0.8974830 0.324 0.100  1.678279e-44       0     Gal
##         mean.in.cluster mean.out.of.cluster
## Pcp4          1.6811115           0.4801027
## Tac1          2.0261848           1.0531066
## Marcks        1.5559093           0.7591471
## Adcyap1       1.0159140           0.3406419
## Nrsn1         1.7072675           0.9505800
## Gal           0.7620386           0.2084508
```

## Finishing up clusters.

At this point in time you should use the tree, markers, domain knowledge, and goals to finalize your clusters. This may mean adjusting PCA to use, mergers clusters together, choosing a new resolutions, etc. When finished you can further name it cluster by something more informative. Ex.

```r
experiment.clusters <- experiment.aggregate
experiment.clusters <- RenameIdents(
  object = experiment.clusters,
  '0' = 'cell_type_A',
  '1' = 'cell_type_B',
  '2' = 'cell_type_C'
)
# and so on

DimPlot(object = experiment.clusters, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

```r
experiment.merged$finalcluster <- Idents(experiment.merged)
```

## Subsetting samples
If you want to look at the representation of just one sample, or sets of samples

```r
experiment.sample2 <- subset(experiment.merged, orig.ident == "UCD_Supp_VitE")

DimPlot(object = experiment.sample2, group.by = "RNA_snn_res.0.5", pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

```r
FeaturePlot(experiment.sample2, features =c('Calca'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

```r
FeaturePlot(experiment.sample2, features =c('Adcyap1'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-25-3.png)<!-- -->

```r
experiment.batch1 <- subset(experiment.merged, batchid == "Batch1")

DimPlot(object = experiment.batch1, group.by = "RNA_snn_res.0.5", pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-25-4.png)<!-- -->

### Adding in a new metadata column representing samples within clusters


```r
experiment.merged$samplecluster = paste(experiment.merged$orig.ident,experiment.merged$finalcluster,sep = '-')

# set the identity to the new variable 
Idents(experiment.merged) <- "samplecluster"

markers.comp <- FindMarkers(experiment.merged, ident.1 = "UCD_Adj_VitE-0", ident.2= c("UCD_Supp_VitE-0","UCD_VitE_Def-0"))

markers.comp
```

```
##                 p_val  avg_logFC pct.1 pct.2    p_val_adj
## Actb     1.183758e-10 -0.5181553 0.993 0.972 1.516513e-06
## Ndufa3   1.726235e-07  0.2771138 0.633 0.419 2.211479e-03
## Pcp4     3.185979e-07  0.3935550 0.745 0.575 4.081557e-03
## Tmsb10   4.479867e-07  0.2762000 0.948 0.897 5.739157e-03
## Rpl21    6.119361e-07  0.3487513 0.861 0.716 7.839513e-03
## Atp5e    1.245447e-05  0.2579366 0.727 0.570 1.595543e-01
## Erp29    5.928481e-05  0.2908362 0.360 0.224 7.594977e-01
## Dbpht2   1.636614e-04  0.2793749 0.303 0.181 1.000000e+00
## Arhgap15 2.074370e-04  0.2892275 0.180 0.089 1.000000e+00
## Gpx3     2.216773e-04  0.2730416 0.382 0.245 1.000000e+00
## Wdfy1    1.831103e-03 -0.2644521 0.056 0.126 1.000000e+00
## Itga6    3.983884e-03  0.2519751 0.112 0.055 1.000000e+00
## Cbx3     4.043714e-03 -0.4100692 0.345 0.407 1.000000e+00
## Cartpt   1.248596e-02  0.2926574 0.172 0.108 1.000000e+00
## Mt2      1.512022e-02  0.2709658 0.187 0.123 1.000000e+00
## Aff3     2.886276e-02 -0.2680419 0.105 0.156 1.000000e+00
## Abhd2    3.420656e-02 -0.4601006 0.341 0.387 1.000000e+00
## Birc6    3.474076e-02 -0.3123984 0.139 0.190 1.000000e+00
## Plp1     4.093181e-02 -0.2754103 0.127 0.181 1.000000e+00
## Hdlbp    4.918260e-02 -0.3744042 0.146 0.188 1.000000e+00
## Cadm3    8.110966e-02 -0.2948876 0.210 0.247 1.000000e+00
## Clasp2   8.530859e-02 -0.3256511 0.195 0.234 1.000000e+00
## Lasp1    9.134266e-02 -0.2663812 0.165 0.204 1.000000e+00
## Usp22    9.683177e-02 -0.2661925 0.288 0.325 1.000000e+00
## Akap9    9.685074e-02 -0.2512918 0.161 0.201 1.000000e+00
## Necab1   1.067273e-01 -0.3338435 0.131 0.167 1.000000e+00
## Zhx1     1.369138e-01 -0.2584958 0.101 0.131 1.000000e+00
## Fam168b  1.679998e-01 -0.2650448 0.176 0.202 1.000000e+00
## Bhlhe41  1.791974e-01 -0.2964435 0.348 0.377 1.000000e+00
## Aqp1     1.877121e-01 -0.2666116 0.356 0.380 1.000000e+00
## Jup      1.933382e-01 -0.2724179 0.315 0.332 1.000000e+00
## Ddhd1    2.063209e-01 -0.2885703 0.109 0.133 1.000000e+00
## Gpsm3    2.149255e-01 -0.2692685 0.097 0.121 1.000000e+00
## Mt1      2.226626e-01  0.2690128 0.464 0.407 1.000000e+00
## Armc8    2.539710e-01 -0.2711362 0.172 0.195 1.000000e+00
## Rtcb     2.887037e-01 -0.2647632 0.154 0.171 1.000000e+00
## Pam      3.206657e-01 -0.3743205 0.509 0.490 1.000000e+00
## Ythdf2   3.773664e-01 -0.2782778 0.311 0.307 1.000000e+00
## Wtap     4.772890e-01 -0.2980818 0.348 0.343 1.000000e+00
## Spry2    5.154505e-01 -0.2640849 0.213 0.217 1.000000e+00
## Acpp     5.297953e-01 -0.2583925 0.195 0.197 1.000000e+00
## Nsg2     6.207279e-01 -0.2622005 0.213 0.211 1.000000e+00
## Setd3    8.494933e-01 -0.2720362 0.285 0.263 1.000000e+00
```

```r
experiment.subset <- subset(experiment.merged, samplecluster %in%  c( "UCD_Adj_VitE-0", "UCD_Supp_VitE-0" ))
DoHeatmap(experiment.subset, features = rownames(markers.comp))
```

```
## Warning in DoHeatmap(experiment.subset, features = rownames(markers.comp)):
## The following features were omitted as they were not found in the
## scale.data slot for the RNA assay: Setd3, Nsg2, Wtap, Ythdf2, Rtcb, Armc8,
## Ddhd1, Jup, Aqp1, Fam168b, Zhx1, Necab1, Akap9, Usp22, Lasp1, Clasp2,
## Hdlbp, Birc6, Aff3, Itga6, Wdfy1, Arhgap15, Erp29, Atp5e, Rpl21, Tmsb10,
## Ndufa3
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

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
