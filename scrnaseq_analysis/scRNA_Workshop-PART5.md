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

experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = use.pcs)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
experiment.aggregate <- FindClusters(
    object = experiment.aggregate, 
    resolution = seq(0.5,4,0.5), 
    verbose = FALSE
)
```

Lets first investigate how many clusters each resolution produces and set it to the smallest resolutions of 0.5 (fewest clusters). 


```r
sapply(grep("res",colnames(experiment.aggregate@meta.data),value = TRUE),
       function(x) length(unique(experiment.aggregate@meta.data[,x])))
```

```
## RNA_snn_res.0.5   RNA_snn_res.1 RNA_snn_res.1.5   RNA_snn_res.2 
##              14              16              20              25 
## RNA_snn_res.2.5   RNA_snn_res.3 RNA_snn_res.3.5   RNA_snn_res.4 
##              26              27              28              29
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
##   1           121           143          106
##   2            77           115          155
##   3            93            82           94
##   4            55            77           86
##   5            51            79           62
##   6            60            66           65
##   7            48            50           45
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

Plot TSNE coloring by the slot 'orig.ident' (sample names).

```r
DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Plot TSNE coloring by the clustering resolution 4

```r
DimPlot(object = experiment.aggregate, group.by="RNA_snn_res.4", pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

FeaturePlot can be used to color cells with a 'feature', non categorical data, like number of UMIs

```r
FeaturePlot(experiment.aggregate, features = c('nCount_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-10-1.png)<!-- -->
and number of genes present

```r
FeaturePlot(experiment.aggregate, features = c('nFeature_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

percent mitochondrial 

```r
FeaturePlot(experiment.aggregate, features = c('percent.mito'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

TSNE plot by cell cycle

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, group.by = "cell.cycle", reduction = "tsne" )
```

<div class="figure" style="text-align: center">
<img src="scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-13-1.png" alt="TSNE Plot by Cell Cycle, No Adjustment"  />
<p class="caption">TSNE Plot by Cell Cycle, No Adjustment</p>
</div>


## Building  a  tree relating the 'average' cell from each cluster. Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space.


```r
experiment.aggregate <- BuildClusterTree(
  experiment.aggregate, dims = use.pcs)

PlotClusterTree(experiment.aggregate)
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


```r
DimPlot(object = experiment.aggregate, pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


```r
experiment.merged <- RenameIdents(
  object = experiment.aggregate,
  '8' = '0'
)
experiment.merged <- RenameIdents(
  object = experiment.merged,
  '6' = '0'
)
DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

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
## Calculating cluster 7
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
## [1] 4950    7
```

```r
head(markers_all)
```

```
##               p_val avg_logFC pct.1 pct.2    p_val_adj cluster   gene
## Nrsn1  6.564989e-95 1.0238837 0.794 0.458 8.410407e-91       0  Nrsn1
## Tac1   2.875256e-69 0.7105631 0.770 0.450 3.683490e-65       0   Tac1
## Stmn1  2.081935e-67 0.6332682 0.868 0.779 2.667167e-63       0  Stmn1
## Marcks 3.076577e-64 0.8034419 0.696 0.397 3.941402e-60       0 Marcks
## Gal    1.020738e-53 0.9448504 0.338 0.099 1.307668e-49       0    Gal
## Calca  2.300792e-53 0.6625024 0.737 0.468 2.947545e-49       0  Calca
```

```r
table(table(markers_all$gene))
```

```
## 
##    1    2    3    4    5 
## 1613  889  390   91    5
```

```r
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

dim(markers_all_single)
```

```
## [1] 1613    7
```

```r
table(table(markers_all_single$gene))
```

```
## 
##    1 
## 1613
```

```r
table(markers_all_single$cluster)
```

```
## 
##   0   1   2   3   4   5   7   9  10  11  12  13 
##  14  86 239 128 145 111 401  65 184  16 136  88
```

```r
head(markers_all_single)
```

```
##                 p_val avg_logFC pct.1 pct.2    p_val_adj cluster     gene
## Nrsn1    6.564989e-95 1.0238837 0.794 0.458 8.410407e-91       0    Nrsn1
## Gm13889  1.355705e-47 0.7209084 0.682 0.484 1.736793e-43       0  Gm13889
## Samsn1   2.815335e-26 0.6310452 0.265 0.114 3.606725e-22       0   Samsn1
## Map1lc3a 1.192513e-21 0.2577715 0.914 0.865 1.527729e-17       0 Map1lc3a
## Nrip1    2.100364e-15 0.6367865 0.316 0.207 2.690776e-11       0    Nrip1
## H2afy    1.675521e-10 0.4126248 0.466 0.384 2.146510e-06       0    H2afy
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
## [1] 60  7
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
## slot for the RNA assay: mt-Nd5, Nwd2, Ctnnd2, Nrip1, Samsn1, Gm13889, Nrsn1
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
##               p_val avg_logFC pct.1 pct.2    p_val_adj cluster   gene
## Nrsn1  6.564989e-95 1.0238837 0.794 0.458 8.410407e-91       0  Nrsn1
## Tac1   2.875256e-69 0.7105631 0.770 0.450 3.683490e-65       0   Tac1
## Stmn1  2.081935e-67 0.6332682 0.868 0.779 2.667167e-63       0  Stmn1
## Marcks 3.076577e-64 0.8034419 0.696 0.397 3.941402e-60       0 Marcks
## Gal    1.020738e-53 0.9448504 0.338 0.099 1.307668e-49       0    Gal
## Calca  2.300792e-53 0.6625024 0.737 0.468 2.947545e-49       0  Calca
##        mean.in.cluster mean.out.of.cluster
## Nrsn1        1.9321788           0.8726252
## Tac1         2.0982494           1.0435836
## Stmn1        2.5283721           1.8028816
## Marcks       1.5974963           0.7586272
## Gal          0.7977859           0.2052266
## Calca        2.4531006           1.3263797
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
##                 p_val  avg_logFC pct.1 pct.2    p_val_adj
## Tmsb10   7.684690e-10  0.2638015 0.988 0.942 9.844857e-06
## Actb     1.321392e-09 -0.5577055 0.996 0.978 1.692836e-05
## Pcp4     8.193620e-09  0.4984027 0.667 0.455 1.049685e-04
## Rpl21    3.063727e-08  0.3658117 0.904 0.755 3.924940e-04
## Rpl23a   4.355414e-07  0.2641398 0.929 0.802 5.579721e-03
## Car8     9.396879e-07  0.2692177 0.129 0.036 1.203834e-02
## Atp5e    1.192436e-06  0.3220314 0.738 0.580 1.527630e-02
## Ndufa3   1.950435e-06  0.2673952 0.617 0.429 2.498703e-02
## Rpl17    5.814455e-06  0.2504313 0.900 0.749 7.448898e-02
## Ap2a1    1.115988e-05  0.2523758 0.254 0.127 1.429693e-01
## Pdap1    1.703453e-05  0.2862285 0.442 0.275 2.182293e-01
## Fxyd7    4.341387e-05  0.3332908 0.458 0.311 5.561751e-01
## Gpx3     4.996700e-05  0.3328014 0.396 0.244 6.401273e-01
## Rhob     7.881333e-05  0.2814958 0.346 0.205 1.000000e+00
## Zfp467   9.421480e-05  0.2834557 0.104 0.035 1.000000e+00
## Mt3      1.395953e-04  0.2602589 0.883 0.749 1.000000e+00
## Dbpht2   1.950345e-04  0.3190892 0.283 0.164 1.000000e+00
## Arhgap15 2.011914e-04  0.3278610 0.188 0.093 1.000000e+00
## Snrpg    2.279383e-04  0.2552111 0.229 0.125 1.000000e+00
## S100a11  4.735212e-04  0.2534807 0.333 0.216 1.000000e+00
## Cd82     1.042477e-03 -0.2966794 0.050 0.125 1.000000e+00
## Erp29    1.091727e-03  0.2758046 0.392 0.276 1.000000e+00
## Gnb1     1.161127e-03 -0.3288667 0.858 0.836 1.000000e+00
## Cbx3     1.520229e-03 -0.4329738 0.321 0.407 1.000000e+00
## BC005624 1.971834e-03  0.2670509 0.267 0.175 1.000000e+00
## mt-Nd4l  2.080357e-03  0.2543612 0.242 0.151 1.000000e+00
## Etv1     2.402617e-03 -0.4891759 0.100 0.180 1.000000e+00
## Cartpt   2.699350e-03  0.3357466 0.196 0.113 1.000000e+00
## Tmem255a 3.556728e-03  0.2540950 0.412 0.305 1.000000e+00
## Mt2      7.173504e-03  0.3353622 0.167 0.098 1.000000e+00
## Birc6    7.707557e-03 -0.3134252 0.133 0.209 1.000000e+00
## Eif5     1.060473e-02 -0.2940629 0.496 0.524 1.000000e+00
## Vdac1    1.336564e-02 -0.2961506 0.488 0.522 1.000000e+00
## Selenop  1.793136e-02 -0.2508145 0.050 0.102 1.000000e+00
## Fam168b  2.078477e-02 -0.2741915 0.112 0.171 1.000000e+00
## Plp1     2.615853e-02 -0.3137040 0.112 0.173 1.000000e+00
## Usp22    3.119112e-02 -0.2925040 0.238 0.298 1.000000e+00
## Aplp2    4.412238e-02 -0.3018389 0.621 0.605 1.000000e+00
## Hdgf     4.610414e-02 -0.2786777 0.308 0.362 1.000000e+00
## Crip1    4.709299e-02 -0.2725004 0.150 0.209 1.000000e+00
## Agap3    4.904363e-02 -0.2512416 0.125 0.176 1.000000e+00
## Hdlbp    5.304203e-02 -0.3205958 0.096 0.140 1.000000e+00
## Phf20    6.476087e-02 -0.2529297 0.146 0.193 1.000000e+00
## Asap1    6.681940e-02 -0.2550599 0.292 0.342 1.000000e+00
## Icmt     6.908485e-02 -0.2824116 0.088 0.127 1.000000e+00
## Necab1   7.668664e-02 -0.3181252 0.146 0.191 1.000000e+00
## Acpp     8.401520e-02 -0.2923132 0.300 0.331 1.000000e+00
## Kcnmb1   1.178421e-01 -0.3334086 0.125 0.160 1.000000e+00
## Nsg2     1.238513e-01 -0.2984264 0.262 0.289 1.000000e+00
## Tle4     1.309881e-01 -0.3080768 0.162 0.191 1.000000e+00
## Mrgpra3  1.318399e-01 -0.3149779 0.092 0.125 1.000000e+00
## Snx7     1.402121e-01 -0.2953552 0.283 0.313 1.000000e+00
## Eno2     1.410445e-01 -0.2647290 0.338 0.360 1.000000e+00
## Abhd2    1.534475e-01 -0.4438863 0.362 0.384 1.000000e+00
## Arf3     1.538958e-01 -0.2711061 0.512 0.505 1.000000e+00
## Mt1      1.725636e-01  0.3288452 0.433 0.378 1.000000e+00
## Clasp2   1.825634e-01 -0.2683631 0.196 0.229 1.000000e+00
## Gpsm3    2.157754e-01 -0.2512666 0.100 0.125 1.000000e+00
## Spry2    2.555499e-01 -0.3280678 0.183 0.205 1.000000e+00
## Ythdf2   2.624758e-01 -0.2975725 0.317 0.320 1.000000e+00
## Bhlhe41  2.718756e-01 -0.2696067 0.379 0.398 1.000000e+00
## Wtap     3.036890e-01 -0.3226680 0.333 0.349 1.000000e+00
## Plxnc1   3.174977e-01 -0.2691579 0.146 0.167 1.000000e+00
## Anp32e   3.453084e-01 -0.2654800 0.183 0.202 1.000000e+00
## Prkca    3.497289e-01 -0.2507596 0.129 0.147 1.000000e+00
## Parp1    4.132646e-01 -0.2515310 0.121 0.136 1.000000e+00
## Setd3    4.297774e-01 -0.2532472 0.246 0.256 1.000000e+00
## Ncam1    4.312967e-01 -0.2596636 0.412 0.398 1.000000e+00
## Pam      4.529740e-01 -0.3829036 0.521 0.505 1.000000e+00
## Ddhd1    5.796372e-01 -0.2632687 0.138 0.144 1.000000e+00
## AI593442 7.259410e-01 -0.2665855 0.200 0.198 1.000000e+00
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
## assay: Ddhd1, Ncam1, Setd3, Anp32e, Wtap, Ythdf2, Clasp2, Arf3, Tle4,
## Nsg2, Necab1, Icmt, Asap1, Phf20, Hdlbp, Hdgf, Usp22, Fam168b, Vdac1, Eif5,
## Birc6, mt-Nd4l, BC005624, Erp29, Snrpg, Arhgap15, Zfp467, Rhob, Pdap1,
## Ap2a1, Rpl17, Ndufa3, Atp5e, Rpl23a, Rpl21, Tmsb10
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
##  [1] nlme_3.1-140        tsne_0.1-3          bitops_1.0-6       
##  [4] RColorBrewer_1.1-2  httr_1.4.0          sctransform_0.2.0  
##  [7] tools_3.6.0         R6_2.4.0            irlba_2.3.3        
## [10] KernSmooth_2.23-15  lazyeval_0.2.2      colorspace_1.4-1   
## [13] npsurv_0.4-0        withr_2.1.2         tidyselect_0.2.5   
## [16] gridExtra_2.3       compiler_3.6.0      plotly_4.9.0       
## [19] labeling_0.3        caTools_1.17.1.2    scales_1.0.0       
## [22] lmtest_0.9-37       ggridges_0.5.1      pbapply_1.4-0      
## [25] stringr_1.4.0       digest_0.6.19       rmarkdown_1.13     
## [28] R.utils_2.9.0       pkgconfig_2.0.2     htmltools_0.3.6    
## [31] bibtex_0.4.2        highr_0.8           htmlwidgets_1.3    
## [34] rlang_0.3.4         zoo_1.8-6           jsonlite_1.6       
## [37] ica_1.0-2           gtools_3.8.1        R.oo_1.22.0        
## [40] magrittr_1.5        Matrix_1.2-17       Rcpp_1.0.1         
## [43] munsell_0.5.0       ape_5.3             reticulate_1.12    
## [46] R.methodsS3_1.7.1   stringi_1.4.3       yaml_2.2.0         
## [49] gbRd_0.4-11         MASS_7.3-51.4       gplots_3.0.1.1     
## [52] Rtsne_0.15          plyr_1.8.4          grid_3.6.0         
## [55] parallel_3.6.0      gdata_2.18.0        listenv_0.7.0      
## [58] ggrepel_0.8.1       crayon_1.3.4        lattice_0.20-38    
## [61] cowplot_0.9.4       splines_3.6.0       SDMTools_1.1-221.1 
## [64] knitr_1.23          pillar_1.4.1        igraph_1.2.4.1     
## [67] future.apply_1.3.0  reshape2_1.4.3      codetools_0.2-16   
## [70] glue_1.3.1          evaluate_0.14       lsei_1.2-0         
## [73] metap_1.1           data.table_1.12.2   png_0.1-7          
## [76] Rdpack_0.11-0       gtable_0.3.0        RANN_2.6.1         
## [79] purrr_0.3.2         tidyr_0.8.3         future_1.13.0      
## [82] assertthat_0.2.1    xfun_0.7            rsvd_1.0.1         
## [85] survival_2.44-1.1   viridisLite_0.3.0   tibble_2.1.3       
## [88] cluster_2.1.0       globals_0.12.4      fitdistrplus_1.0-14
## [91] ROCR_1.0-7
```
