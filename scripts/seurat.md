R Notebook
================

### Task specification: Case study

In the following tutorial, a dataset of Peripheral Blood Mononuclear
Cells(PBMC) is analyzed. The dataset is from 10X Genomics and consists
out of 2700 single cells sequenced with Illumina’s NextSeq 500.

### Install required packages

install.packages(‘Seurat’) install.packages(‘umap’)

### Loading libraries

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.0.5

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(Seurat)
```

    ## Warning: package 'Seurat' was built under R version 4.0.5

    ## Attaching SeuratObject

``` r
library(patchwork)
```

    ## Warning: package 'patchwork' was built under R version 4.0.3

``` r
# library(umap)
```

## Setup the seurat Object

### Reading in the data: The Read10X() function is used to read the output of the Cellranger pipeline from 10X. Cellranger is an analysis pipeline that processes single-cell data to align reads, generates feature barcode matrices, performs clustering and other secondary analysis. The output is an UMI(Unique Molecular Identifier) count martix. UMIs are short sequences that are added to DNA sequences to uniquely identify the DNA molecule. The values in the matrix represent the number of molecules for each feature (gene, row) detected in each cell (column). The Count Matrix is used to create the Seurat object, which contains both the data and the analyses, e.g. PCA, clustering results, for a data set.

``` r
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
pbmc
```

    ## An object of class Seurat 
    ## 13714 features across 2700 samples within 1 assay 
    ## Active assay: RNA (13714 features, 0 variable features)

## Preprocessing and Quality control

### Preprocessing means the selection and filtering of cells based on QC metrices, normalization of the data and scaling or the detection of highly variable features. QC metrices are defined by the user and are, for example: the number of unique genes in each cell (cells with poor quality or empty droplets often contain only a few genes, duplicates or multiplets can have very high gene counts). Unique genes correlate strongly with the number of molecules that are detected within a cell. Another metric is the percentage of reads that map to the mitochondrial genome. Dying and poor quality cells often have extensive mitochondrial contamination. The calculation is done with PercentageFeatureSet (), a function that calculates the percentage of counts from a certain feature set. (all genes starting with MT).

``` r
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

## Visualization of QC Metrices

### The filter criteria for the cells are defined as follows: Unique feature count &gt; 2500 or &lt; 200 and mitochondrial counts &gt; 5%.

``` r
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

![](seurat_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> \#\#
Visualization of feature-feature relationships by FeatureScatter

``` r
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

![](seurat_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Subsets with specific filter criteria can also be defined

``` r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

## Normalization of the data

### After the filter steps, the data must be normalized. By default, the “LogNormalize” method is used for this, which normalizes the feature expression measurements for each cell by the total expression, multiplies by the scale factor (default 10000) and log-transforms the result.

``` r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

``` r
### the same
pbmc <- NormalizeData(pbmc)
```

## Feature selection for the identification of highly variable features

### A subset of features is defined with high cell-to-cell variation in the expression (highly expressed in some cells and lowly expressed in others); the further focus of the analysis is then placed on these genes. Therefore the FindVariablesFeatures () function is used which allows directy modeling of the mean-variance relationship in single-cell data. Default: 2000 features per dataset, used e.g. for PCA

``` r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results

``` r
plot1 + plot2
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](seurat_files/figure-gfm/unnamed-chunk-9-1.png)<!-- --> \#\# Scaling
of the data \#\#\# Scaling is the linear transformation which is the
standard before dimension reduction techniques like PCA. (mean
expression across cells is 0 and variance across cells is 1). Scaling is
intended to achieve equal weighting in downstream analysis so that
highly expressed genes do not dominate.

``` r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

    ## Centering and scaling data matrix

## Linear dimension reduction

### The next step is to perform PCA on the scaled data. By default, the previously defined features are used, but a new subset can also be defined.

``` r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

    ## PC_ 1 
    ## Positive:  CST3, TYROBP, LST1, AIF1, FTL, FTH1, LYZ, FCN1, S100A9, TYMP 
    ##     FCER1G, CFD, LGALS1, S100A8, CTSS, LGALS2, SERPINA1, IFITM3, SPI1, CFP 
    ##     PSAP, IFI30, SAT1, COTL1, S100A11, NPC2, GRN, LGALS3, GSTP1, PYCARD 
    ## Negative:  MALAT1, LTB, IL32, IL7R, CD2, B2M, ACAP1, CD27, STK17A, CTSW 
    ##     CD247, GIMAP5, AQP3, CCL5, SELL, TRAF3IP3, GZMA, MAL, CST7, ITM2A 
    ##     MYC, GIMAP7, HOPX, BEX2, LDLRAP1, GZMK, ETS1, ZAP70, TNFAIP8, RIC3 
    ## PC_ 2 
    ## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1, HLA-DRA, LINC00926, CD79B, HLA-DRB1, CD74 
    ##     HLA-DMA, HLA-DPB1, HLA-DQA2, CD37, HLA-DRB5, HLA-DMB, HLA-DPA1, FCRLA, HVCN1, LTB 
    ##     BLNK, P2RX5, IGLL5, IRF8, SWAP70, ARHGAP24, FCGR2B, SMIM14, PPP1R14A, C16orf74 
    ## Negative:  NKG7, PRF1, CST7, GZMB, GZMA, FGFBP2, CTSW, GNLY, B2M, SPON2 
    ##     CCL4, GZMH, FCGR3A, CCL5, CD247, XCL2, CLIC3, AKR1C3, SRGN, HOPX 
    ##     TTC38, APMAP, CTSC, S100A4, IGFBP7, ANXA1, ID2, IL32, XCL1, RHOC 
    ## PC_ 3 
    ## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1, HLA-DPA1, CD74, MS4A1, HLA-DRB1, HLA-DRA 
    ##     HLA-DRB5, HLA-DQA2, TCL1A, LINC00926, HLA-DMB, HLA-DMA, CD37, HVCN1, FCRLA, IRF8 
    ##     PLAC8, BLNK, MALAT1, SMIM14, PLD4, LAT2, IGLL5, P2RX5, SWAP70, FCGR2B 
    ## Negative:  PPBP, PF4, SDPR, SPARC, GNG11, NRGN, GP9, RGS18, TUBB1, CLU 
    ##     HIST1H2AC, AP001189.4, ITGA2B, CD9, TMEM40, PTCRA, CA2, ACRBP, MMD, TREML1 
    ##     NGFRAP1, F13A1, SEPT5, RUFY1, TSC22D1, MPP1, CMTM5, RP11-367G6.3, MYL9, GP1BA 
    ## PC_ 4 
    ## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1, CD74, HLA-DPB1, HIST1H2AC, PF4, TCL1A 
    ##     SDPR, HLA-DPA1, HLA-DRB1, HLA-DQA2, HLA-DRA, PPBP, LINC00926, GNG11, HLA-DRB5, SPARC 
    ##     GP9, AP001189.4, CA2, PTCRA, CD9, NRGN, RGS18, GZMB, CLU, TUBB1 
    ## Negative:  VIM, IL7R, S100A6, IL32, S100A8, S100A4, GIMAP7, S100A10, S100A9, MAL 
    ##     AQP3, CD2, CD14, FYB, LGALS2, GIMAP4, ANXA1, CD27, FCN1, RBP7 
    ##     LYZ, S100A11, GIMAP5, MS4A6A, S100A12, FOLR3, TRABD2A, AIF1, IL8, IFI6 
    ## PC_ 5 
    ## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY, CCL4, CST7, PRF1, GZMA, SPON2 
    ##     GZMH, S100A9, LGALS2, CCL3, CTSW, XCL2, CD14, CLIC3, S100A12, CCL5 
    ##     RBP7, MS4A6A, GSTP1, FOLR3, IGFBP7, TYROBP, TTC38, AKR1C3, XCL1, HOPX 
    ## Negative:  LTB, IL7R, CKB, VIM, MS4A7, AQP3, CYTIP, RP11-290F20.3, SIGLEC10, HMOX1 
    ##     PTGES3, LILRB2, MAL, CD27, HN1, CD2, GDI2, ANXA5, CORO1B, TUBA1B 
    ##     FAM110A, ATP1A1, TRADD, PPA1, CCDC109B, ABRACL, CTD-2006K23.1, WARS, VMO1, FYB

## Visualizing cells and features that define PCA

``` r
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
```

    ## PC_ 1 
    ## Positive:  CST3, TYROBP, LST1, AIF1, FTL 
    ## Negative:  MALAT1, LTB, IL32, IL7R, CD2 
    ## PC_ 2 
    ## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1 
    ## Negative:  NKG7, PRF1, CST7, GZMB, GZMA 
    ## PC_ 3 
    ## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1 
    ## Negative:  PPBP, PF4, SDPR, SPARC, GNG11 
    ## PC_ 4 
    ## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1 
    ## Negative:  VIM, IL7R, S100A6, IL32, S100A8 
    ## PC_ 5 
    ## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY 
    ## Negative:  LTB, IL7R, CKB, VIM, MS4A7

``` r
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```

![](seurat_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
DimPlot(pbmc, reduction = "pca")
```

![](seurat_files/figure-gfm/unnamed-chunk-14-1.png)<!-- --> \#\#
Visualization with DimHeatmap() \#\#\# DimHeatmap() allows visualization
of the primary source of heterogeneity in the dataset and can help to
decide which PCs to use for further downstream analysis. Cells and
features are ordered by their PCA scores, which also increases speed.
this is especially important for large datasets.

``` r
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
```

![](seurat_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```

![](seurat_files/figure-gfm/unnamed-chunk-16-1.png)<!-- --> \#\#
Determine the dimensionality of a dataset \#\#\# To eliminate technical
noise in the individual features, Seurat clusters the cells based on
their PCA scores, each PC representing a so-called “metafeature” which
combines information about a correlated feature set. In a resampling
test, a subset of the data (1% by default) is randomly selected and the
PCA is repeated. A null distribution of feature scores is identified…the
procedure is repeated…This allows one to identify significant PCs that
contain many low p-values.

``` r
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
```

## Visualization

### Comparison of the distribution of p-values for each PC with a unit distribution. Significant PCs show a strong enrichment of features with low p-values (solid curve above the dashed line), we can see a sharp drop-off in significance after the first 10 - 12 PCs.

``` r
JackStrawPlot(pbmc, dims = 1:15)
```

    ## Warning: Removed 23504 rows containing missing values (geom_point).

![](seurat_files/figure-gfm/unnamed-chunk-18-1.png)<!-- --> \#\#
Alternative: Elbow plot \#\#\# This plot shows a ranking of PCs based on
the percentage of variance explained. Cutoff beteen PC 7 -12

``` r
ElbowPlot(pbmc)
```

![](seurat_files/figure-gfm/unnamed-chunk-19-1.png)<!-- --> \#\#\#
Conclusio: To identify the true dimensionality of a dataset is a
challenge. Therefore, several approaches should always be included. 1.
Supervised (exploring PCs to determine sources for heterogeneity),
e.g. in combination with GSEA 2. Statistical test based on the random
null model (time consuming for large datasets, no clear PC cutoff) 3.
Heuristic (calculated instantly)

We choose 10 PCS for the cutoff, but the following still needs to be
considered: Dendritic cell and NK aficionados may recognize that genes
strongly associated with PCs 12 and 13 define rare immune subsets
(i.e. MZB1 is a marker for plasmacytoid DCs). However, these groups are
so rare, they are difficult to distinguish from background noise for a
dataset of this size without prior knowledge. &gt; repeat downstream
analyzes with a different number of PCs (10, 15, or even 50!)

## Cluster the cells

### Seutrat uses a graph-based clustering approach. The distance metric remains the same and is based on the previously identified PCs. The approach is based on manuscripts applied to graph-based clustering approaches with scRNA seq data and is characterized by a graph structure, e.g. KNN nearest neighbors. (edges drawn between cells with similar feature expression patterns) and then graph is partitioned into highly interconnected communities.

First, a KNN graph is constructed based on the Euclidean distance (in
PCA space), edge weights are refined between all two cells based on the
overlap in their neighborhood. This step is done using the
FindNeighbors() function with the first 10 PCs as input. To cluster the
cells, modularity optimization techniques such as the Louvain algorithm
or SLM are applied. Groups are iteratively clustered with the goal of
optimizing the standard modularity function. The FindClusters() function
is used for this, this sets the granularity for downstream clustering
(higher values, more clusters). Settings: 0.4 - 1.2 for single-cell
datasets with about 3K cells. The larger the dataset, the more likely is
an optimal solution. With the Idents() function the clusters can be
found.

``` r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2638
    ## Number of edges: 95965
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8723
    ## Number of communities: 9
    ## Elapsed time: 0 seconds

``` r
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```

    ## AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 
    ##                2                3                2                1 
    ## AAACCGTGTATGCG-1 
    ##                6 
    ## Levels: 0 1 2 3 4 5 6 7 8

## Non-linear dimensional reduction(UMAP/tSNE)

### The goal of nonlinear dimensional reduction techniques is to learn underlying manifold of the data and place similar cells together in low-dimensional space. Cells in the previously determined clusters are to be displayed co-localized in the plots. The same PCs are used as input as for the cluster analysis.

``` r
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 20:33:35 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 20:33:35 Read 2638 rows and found 10 numeric columns

    ## 20:33:35 Using Annoy for neighbor search, n_neighbors = 30

    ## 20:33:35 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 20:33:35 Writing NN index file to temp file /tmp/RtmpxBzELd/file6cd064223db7
    ## 20:33:35 Searching Annoy index using 1 thread, search_k = 3000
    ## 20:33:36 Annoy recall = 100%
    ## 20:33:36 Commencing smooth kNN distance calibration using 1 thread
    ## 20:33:37 Initializing from normalized Laplacian + noise
    ## 20:33:37 Commencing optimization for 500 epochs, with 105124 positive edges
    ## 20:33:39 Optimization finished

``` r
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
```

![](seurat_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
### output ordner angelegt
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
```

## Search for differentially expressed traits (cluster biomarkers)

### Now, markers are to be found which define clusters through differential gene expression. By default, positive and negative markers of a single cluster are identified compared to all other cells. Set threshold value: A trait must be recognized to a minimum percentage in one of the two cell groups and must be expressed differently between the two groups. (both set to 0 leads to dramatic time savings). The maximum number of cells can also be defined (increases the speed)

``` r
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
```

    ## For a more efficient implementation of the Wilcoxon Rank Sum Test,
    ## (default method for FindMarkers) please install the limma package
    ## --------------------------------------------
    ## install.packages('BiocManager')
    ## BiocManager::install('limma')
    ## --------------------------------------------
    ## After installation of limma, Seurat will automatically use the more 
    ## efficient implementation (no further action necessary).
    ## This message will be shown once per session

``` r
head(cluster2.markers, n = 5)
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## IL32 2.593535e-91  1.2154360 0.949 0.466 3.556774e-87
    ## LTB  7.994465e-87  1.2828597 0.981 0.644 1.096361e-82
    ## CD3D 3.922451e-70  0.9359210 0.922 0.433 5.379250e-66
    ## IL7R 1.130870e-66  1.1776027 0.748 0.327 1.550876e-62
    ## LDHB 4.082189e-65  0.8837324 0.953 0.614 5.598314e-61

``` r
### find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
```

    ##                       p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## FCGR3A        2.150929e-209   4.267579 0.975 0.039 2.949784e-205
    ## IFITM3        6.103366e-199   3.877105 0.975 0.048 8.370156e-195
    ## CFD           8.891428e-198   3.411039 0.938 0.037 1.219370e-193
    ## CD68          2.374425e-194   3.014535 0.926 0.035 3.256286e-190
    ## RP11-290F20.3 9.308287e-191   2.722684 0.840 0.016 1.276538e-186

``` r
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

    ## Calculating cluster 0

    ## Calculating cluster 1

    ## Calculating cluster 2

    ## Calculating cluster 3

    ## Calculating cluster 4

    ## Calculating cluster 5

    ## Calculating cluster 6

    ## Calculating cluster 7

    ## Calculating cluster 8

``` r
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
```

    ## Registered S3 method overwritten by 'cli':
    ##   method     from         
    ##   print.boxx spatstat.geom

    ## # A tibble: 18 × 7
    ## # Groups:   cluster [9]
    ##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    
    ##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
    ##  1 1.17e- 83       1.33 0.435 0.108 1.60e- 79 0       CCR7    
    ##  2 1.74e-109       1.07 0.897 0.593 2.39e-105 0       LDHB    
    ##  3 0               5.57 0.996 0.215 0         1       S100A9  
    ##  4 0               5.48 0.975 0.121 0         1       S100A8  
    ##  5 7.99e- 87       1.28 0.981 0.644 1.10e- 82 2       LTB     
    ##  6 2.61e- 59       1.24 0.424 0.111 3.58e- 55 2       AQP3    
    ##  7 0               4.31 0.936 0.041 0         3       CD79A   
    ##  8 9.48e-271       3.59 0.622 0.022 1.30e-266 3       TCL1A   
    ##  9 4.93e-169       3.01 0.595 0.056 6.76e-165 4       GZMK    
    ## 10 1.17e-178       2.97 0.957 0.241 1.60e-174 4       CCL5    
    ## 11 3.51e-184       3.31 0.975 0.134 4.82e-180 5       FCGR3A  
    ## 12 2.03e-125       3.09 1     0.315 2.78e-121 5       LST1    
    ## 13 6.82e-175       4.92 0.958 0.135 9.36e-171 6       GNLY    
    ## 14 1.05e-265       4.89 0.986 0.071 1.44e-261 6       GZMB    
    ## 15 1.48e-220       3.87 0.812 0.011 2.03e-216 7       FCER1A  
    ## 16 1.67e- 21       2.87 1     0.513 2.28e- 17 7       HLA-DPB1
    ## 17 3.68e-110       8.58 1     0.024 5.05e-106 8       PPBP    
    ## 18 7.73e-200       7.24 1     0.01  1.06e-195 8       PF4

## Visualization of the marker expression

### The violin plot shows the expression probability distribution over the clusters. The feature plot visualizes the feature expression in a tSNE or PCA plot. RidgePlot (), CellScatter () and DotPlot () can be further options.

``` r
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
```

![](seurat_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
```

![](seurat_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
```

![](seurat_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## DoHeatmap () can be used to generate an expression heatmap for cells and features. The top 20 markers for each cluster are plotted here.

``` r
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

![](seurat_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

## Assigning cell type identity to clusters

### As the last step in the workflow, canonical markers are sent to assign the unbiased clustering to known cell types (9 assignments: “Naive CD4 T”, “CD14 + Mono”, “Memory CD4 T”, “B”, “CD8 T”, " FCGR3A + Mono “,”NK“,”DC“,”Platelet"). The new clusters are visualized as umap with DimPlot ().

``` r
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](seurat_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
```
