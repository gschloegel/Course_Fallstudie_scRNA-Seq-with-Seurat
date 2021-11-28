# Case Study - Replication of Seurat Workflow

## Introduction

replication crisis

questions point 1 scRNA-Seq
advantages

steps

common problems

major challenges

what problems the seurat workflow addresses

## Methods

set up environment
anaconda create env command
exact version of the packages

### work flow description (basic details in R Notebook)

* pre processing
  * select cells, discard low quality cells (very few genes, too many genes - doublets or triplets, high amount of mitochondrial genes)
* Normalize the data
  * "LogNormalize"
* Identification of highly variable features
  * describe method (vst)
  * use only variable features for further analysis
* PCA
  * scale data (mean = 0, variance = 1)
  * PCA on variable features
  * visualize results in different ways (Loadings, scatter plot PCs, heat maps)
  * Heat maps top features for each PC, top 500 cells (color scale is missing)
  * decide with PC(s) to use for further analysis
* determine 'dimensionality' of dataset (PCs used for further analysis)
  * Jack Strow Plot (describe)
  * Elbow plot
* Cluster the cells
  * graph based clustering
  * just use identified PCs
* UMAP: non-linear dimension reduction / display clusters
* finding differentially expressed features
  * use VlnPlot to look at genes with big loadings. See in which clusters these genes are expressed.
  + Feature Plot to see where the features are expressed
* Assigning cell type identity to clusters
  * use markers to determine cell type and label the clusters

mention: compare with rules from Sandve 2013

do analysis on own dataset

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

## Discussion

