# Case Study - Replication of Seurat Workflow

## Introduction

Replication of results is an important part of science. As it is important to avoid the replication crisis in data science, we show here how the Seurat tutorial for single cell RNA sequencing (scRNA-Seq) can be replicated and how it can be adapted to a different dataset.

Single cell RNA Seq is characterized by the fact that each sequencing library represents a single cell. Bulk RNA-Seq used to be standard for all -omics disciplines and represents a whole cell population. Within a single tissue, there are many similar cells that may have specific functions, which can be addressed with scRNA Seq. Bulk sequencing studies the genome/transcriptome, e.g., the differential gene expression of healthy and diseased individuals is compared as an average representation (studying biomarkers and biology of diseases). scRNA Seq has several advantages over bulk RNA-seq, including more accurate prediction of disease, discovery of new cell types, and the ability to study gene expression in individual cells (e.g., finding differences in specific cell types).

scRNA Seq enables detailed analysis of individual cells and reveals cellular heterogeneity. In comparison, the bulk method measures the average expression across a population of cells and the cellular heterogeneity is masked. ScRNA Seq provides high resolution and transcriptional profiling of thousands of individual cells, allows to understand gene expression at a single cell level and finds differences within a heterogenous sample. It measures the distribution of the expression level for each gene across a population of cells. In comparison, a bulk sample represents many cells where all kinds of RNA sequences from the sample are extracted without filtering and enrichment, measuring the average expression level for each gene across a large population of input cells (e.g., a sample from the same tissue from different species).

basic steps - described in methods part

problems and how they are solved

major chalanges and how the are solved

questions point 1 scRNA-Seq
advantages
n
steps

common problems

major challenges

In this text we use the Seurat workflow to robustly separate different cell types in the sample. The target is to integrate the whole workflow from data pre-processing to separation of cell types and the assignment of cell types to the clusters. It used highly variable features to get reliable results that can be compared between different samples.

## Methods

### Configuration steps

First we set up the environment. While there are different options for this we decided to use anaconda (miniconda) for this task. We set up the conda forge channel as default and run the following command

`mamba create -n seurat r-base=4.1.0 r-ggplot2=3.3.5 patchwork=1.1.1 SeuratObject=4.0.2 Seurat=4.0.4 dplyr=1.0.7 r-knitr=1.33`

to create the environment. We kept the package versions used in the tutorial to avoid the risk of incompatibilities and tried to keep the environment as small as possible. The environment was exported as .yaml file to allow easy replication.

mamba is used instead of conda to to speed up the environment management. This is especially relevant as we use the large conda forge channel.

It is also required to create the folder structure manually. The notebooks have to be in a subfolder (scripts in our case), the exact location of the data files has to be set in the notebook and an output folder must be present before running the code.

## Work flow description (basic details in R Notebook)

The workflow runs in the following basic steps

### Pre-processing
First we discard low quality cells, e.g. cells with very few genes, or suspected doublets and triplets with an high amount of genes. We normalize the data using a logarithmic method. 

### Identification of highly variable features and PCA

Next we limit the further steps to the features showing a highly variable. This ensure that we get less noise form the other features. This makes the later calculation faster and more importantly leads to better results.

For the PCA the data from the variable features is normalized (mean = 0, variance = 1). After the PCA is performed the results are visualized. The loadings are presented, a scatter plot of the first PCs is created and head maps for the first PCs are created. This helps to identify where the relevant information can be found and which PCs and which loadings/genes are used later.

#### determine 'dimensionality' of dataset

After the PCA we need to know how many dimensions are needed for the further analyses. In the Tutorial 2 heuristic methods are used. The first is the Jack Straw plot and the second the elbow plot.
  * describe how the plots are used

#### Cluster the cells

A graph based clustering approach is used for the identified relevant PCs. As we do not know a priory how many clusters/cell type we will find this is a sensible approach.

The clusters are displayed using the UMAP non linear dimension reduction.

#### finding differentially expressed features

Finally we want to identify the clusters. We use VlnPlot to see in which cluster the genes with the highest loadings are expressed. Additionally the genes with the highest fold change for each cluster are calculated. With additional information about the cell types expected in the sample the clusters can be assigned to cell types.

### Assess reproducibility of code

cite Sandve 2013

The rules form the Article are used to assess how easy the code can be reproduced and reused.

### Analysis on own dataset

cite source of dataset

A new dataset was introduced to test the possibility to reuse the code.

## Results

### discuss rules from Sandve 2013

1. For Every Result, Keep Track of How It Was Produced

The Tutorial is available as a R-notebook and vignette. So all Steps and R commands are available 

2. Avoid Manual Data Manipulation Steps

The data are not manually manipulated, but there are a couple of manual interventions in the script, for example the cut off in the pre-processing, the number of PCs used for further analyses and the feature selection for the cell type assignment. Without concrete rules for these decisions reproducible results are difficult to obtain.

3. Rule 3: Archive the Exact Versions of All External Programs Used

This is done with the `sessionInfo()` command. The version of R and the used packages is available.

4. Version Control All Custom Scripts

The Seurat package has its own GitHub repository. Past versions of the Tutorial are available.

5. Record All Intermediate Results, When Possible in Standardized Formats

The created Seurat Object is saved as an .rds file. As this is standard in R and further analysis will be done with R as well, this is a good choice.

6. For Analyses That Include Randomness, Note Underlying Random Seeds
results of own dataset (including plots) - what changes

The way the Seurat packages deals with randomness, was confusing for us. The seed is not set in the notebook, but hidden in the Seurat package. This leeds to the strange situation that for example the Jack Straw Plot differs visibly, if not significantly, from the tutorial, but there is no random change when run again as expected. The difference is due to the different handling of seed values on different operating systems.

For us it had been nicer to set the seed in the final notebook to avoid confusion.

7. Always Store Raw Data behind Plots

Most of the plots just require the Seurat object, which is covered in Rule 5. It would be nice to store additional data, like used features and cluster labels' in an easier to read form. In the current form the are in the code. So everything is reproducible, but not always easy to find.

8. Generate Hierarchical Analysis Output, Allowing Layers of Increasing Detail to Be Inspected

The Tutorial is more a way to show the capabilities of the package than tha good way to represent results. Therefore rule is not applied.

9. Connect Textual Statements to Underlying Results

Similar to rule 8 this is missing as a thorough explanation of the results is not part of the tutorial.

10. Provide Public Access to Scripts, Runs, and Results

The code is available in the GitHub repository.

### Discuss own results

summarize what is written in the skript

## Discussion

It is possible to reproduce the tutorial and use the method on a different dataset. While it worked well there are some options to make it easier. First the form of a notebook makes it difficult to see where parameters are set and manual interventions are necessary. It is hard to make changes without checking the whole code. In our opinion it would be easier to write the code in functions and use all variables as function parameters.

It would be possible to do the analysis in 2 well defined steps. Step 1 starting from pre-processing and ending with the analysis of the PCs. At this point manual intervention is done. The rest of the analyses can be performed without manual intervention in a second function or method.