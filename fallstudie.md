# Case Study - Replication of Seurat Workflow

Autoren: Guido Schlögel 00727019
         Sonja Tockner 00708717

## Introduction

Replication of results is an important part of science. As it is important to avoid the replication crisis in data science, we show here how the Seurat tutorial for single cell RNA sequencing (scRNA-Seq) can be replicated and how it can be adapted to a different dataset.

## 1.	Preparation. As a preparation for the data analysis workflow, thoroughly read the given manuscripts (s. literature)

###  a.	Advantages of the scRNA-Seq method compared to bulk RNA-Seq

Single cell RNA Seq is characterized by the fact that each sequencing library represents a single cell. Bulk RNA-Seq used to be standard for all -omics disciplines and represents a whole cell population. Within a single tissue, there are many similar cells that may have specific functions, which can be addressed with scRNA Seq. Bulk sequencing studies the genome/transcriptome, e.g., the differential gene expression of healthy and diseased individuals is compared as an average representation (studying biomarkers and biology of diseases). scRNA Seq has several advantages over bulk RNA-seq, including more accurate prediction of disease, discovery of new cell types, and the ability to study gene expression in individual cells (e.g., finding differences in specific cell types).

scRNA Seq enables detailed analysis of individual cells and reveals cellular heterogeneity. In comparison, the bulk method measures the average expression across a population of cells and the cellular heterogeneity is masked. ScRNA Seq provides high resolution and transcriptional profiling of thousands of individual cells, allows to understand gene expression at a single cell level and finds differences within a heterogenous sample. It measures the distribution of the expression level for each gene across a population of cells. In comparison, a bulk sample represents many cells where all kinds of RNA sequences from the sample are extracted without filtering and enrichment, measuring the average expression level for each gene across a large population of input cells (e.g., a sample from the same tissue from different species).

### b.	Basic steps of the scRNA-Seq analysis

The workflow for scRNA-Seq analysis basically consists of the preprocessing step of the raw data (quality control, normalization, data correction, feature selection and dimensions reduction) and the cell-and gene level downstream analysis of the count data.  
Biological tissue samples are used as input material. First, single-cell dissociation is performed, a suspension is generated, and the tissue is digested. This is followed by the isolation of the single cells. A distinction is made between plate-based and droplet-based methods. The problem can arise that multiple cells are captured together (multiplets, doublets), nonviable cells are captured, or no cells are captured (empty wells).  In the next step, library construction is performed (breaking down the cell membrane), intracellular mRNA is captured, reverse-transcribed to cDNA molecules and amplified. The mRNA is then labeled with barcodes and possibly UMIs. Libraries are pooled for sequencing. Sequencing produces read data and is submitted to quality control. 

•	Preprocessing and Visualization

In this step, matrices of molecular counts or read counts are generated from the raw data. Pipelines such as CellRanger or ZUMIs are used here. Read and count matrix differ in level of measuring noise. The matrix has the dimension of the barcodes multiplied by the number of transcripts. An error that can occur here is that a barcode labels several or even no cells.

Quality control ensures that all cellular barcode data correspond to viable cells. 3 different QC covariates can be considered:
1. Number of counts per barcode (count depth)
2. Number of genes via barcode
3. Fraction of counts from mitochondrial genes per barcode

The results are examined for outliers in the peaks, the filtering is carried out with a set threshold. Cells with unexpectedly high counts and many detected genes may represent doublets. Considering these 3 covariates in isolation can lead to misinterpretation. QC should be considered jointly when univariate threshold decisions are made. The threshold should be set as permissive as possible to filtering out viable cells. QC steps must also be carried out at the transcript level. Raw count matrices often include more than 20000 cells. Filtering eliminates genes expressed only in a few genes that are not informative of the cellular heterogeneity. The choice of the threshold should scale with the number of cells in the dataset. The further quality control takes place directly on the count data. Ambient gene expression means that counts that are not from barcoded cells but from other lysed cells indicate mRNA contamination prior to library construction and must also be corrected. In data sets with low quality, it may be necessary to set stricter thresholds. Furthermore, the QC may have to be readjusted later.
Each count represents a successful capture, reverse transcription, and sequencing of a molecule of the cellular mRNA.

•	Normalization

Normalization is the scaling of the count data to obtain correct relative gene expression abundances between cells. One scRNA-specific normalization method is CPM (counts per million, count depth scaling). CPM assumes that all cells in the dataset initially contained an equal number of mRNA molecules and count depth differences arise only due to sampling. Single cell datasets consist of heterogenous cell populations with varying size and molecule counts, that means that more complex normalization strategies are necessary. An extension to CPM would be to exclude gene levels that account for at least 5% of the total counts in any cell when calculating the size factors. This allows count variability in few higher expressed genes. Scran is another possibility, it is a pooling based size factor estimation method and a top performing normalization method. These 3 use linear, global scale to normalize count data.
Nonlinear methods are for example parametric modeling of the data. These methods can correct technical and biological data (e.g., batch correction for cell cycle effects). In plate-based methods, batch effects between plates must be corrected. ScRNA techniques are divided into full-length (e.g., TPM normalization) and 3' enrichment. Scaling of gene counts improves comparison between genes (zero mean, unit variance) and all genes are weighted equally for downstream analysis. 
After normalization, the data matrices must be log (x+1) normalized to measure changes in expressions. The distance between log-transformed expression values represents the logFold Change. 

•	Data correction and integration

The normalized data may still contain unwanted variability. Therefore, the data must be corrected (for batch, dropout, and cell cycle effects). Count depth and batch are part of the technical effects. Batch effects can occur when cells are traded in distinct groups (e.g., cells on different chips). Batch correction of bulk RNA-Seq corrects the batch effect between samples or cells in the same experiment. The correction is made using linear methods. Data integration (e.g., CCA) means the integration of data from multiple experiments with non-linear methods. The best option is pre-emitting the effect and avoiding it through good experimental design (e.g., pooling cells).
Expression recovery means reducing the noise in the data set, for example by removing dropouts. Data integration and batch corrections should be performed by different methods because data integration tools may overcorrect simple batch effects. 

•	Feature selection, dimensions reduction and visualization

Feature selection means that the data set is filtered, for example by removing zero counts. Only genes that are "informative" of the variability in the data are retained. For the downstream analysis about 1000 to 5000 HVGs (High Variable Genes) are selected. 
Dimension reduction: The visualization describes the dataset in 2 to 3 dimensions, summarization does not prescribe the number of output compounds. A distinction is made between linear and non-linear methods. PCA, a linear method, is the basis for many clustering tools as a preprocessing step. Diffusion maps belong to the nonlinear methods and are used for differentiation, each component highlights the heterogeneity of a different cell population. 
Non-linear methods like tSNE or UMAP are used for visualization. UMAP can scale many cells and can summarize data in more than two dimensions. Reduced data are rather used for summarizing the data, because the differential gene expression cannot be tested properly, for the expression profile the corrected data are used. A problem with data correction is always that the data are over/under corrected.

•	Downstream Analysis

The downstream analysis enables to extract biological insights and describe the underlying biological system by fitting models to the data, e.g., groups of cells with similar gene expression profiles represent cell-type clusters. A distinction is made between cell- and gene-level approaches. The first is clustering, categorization in groups and trajectories. Cells are grouped by similar gene expression profiles (calculation with distance metrices, e. g. Euclidean distance, cosine similarity, Louvain algorithm). Clustering usually takes place with the PCA reduced data and is an unsupervised machine learning method. Cells are divided into K clusters by determining the cluster centroids, and the cells are assigned to the closest cluster centroids.
For cluster annotation clustered data are analyzed by finding the gene signature of each cluster (“marker genes” for biological labeling). Most of the time, the focus is on genes that are upregulated.
The trajectory analysis creates dynamic models for gene expression. Biological processes are continuous, the aim is finding paths through cellular space that minimize transcriptional changes between neighboring cells.
Gene expression dynamics analyze trajectory on gene level. Regulator genes are important for understanding how biological processes are triggered. To represent the static and dynamic nature of the data one can e.g., use the tool PAGA, which represents single cell clusters as nodes and trajectories between cluster as edges.
Gene level analysis examines the molecular signals in the data, a distinction is made between 3 methods:

1. Differential gene expression analysis
This method comes from bulk gene expression analysis and examines the differential gene expression between two different conditions. DESeq2 and EdgeR are tools preferably used, weight estimation (introduce gene weights) can be included. It should be applied to measured data.

2. Gene set analysis
Genes are grouped based on their involvement in common biological processes, e.g., use of paired gene labels to perform ligand-receptor analysis.

3. Gene regulatory networks
Genes do not function independently. The expression level of a gene is determined by a complex interplay of regulatory interactions. Networks make it possible to discover the underlying regulatory interactions and to measure gene co-expression.

### c.	Common problems and how are they typically solved

Actually, the main task is to better understand heterogeneous cell populations. With classical methods it is not possible to get results with good resolution. Therefore, we nees methods that take this problem into account. (e.g. one can better understand tissues, use in cancer research, analysis of microbial systems...). Predefined workflows using specialized software simplify the reproducibility of complex issues.

### d.	Major challenges in integrating single-cell transcriptomic data across different
### conditions, technologies, and species and how they can be solved

ScRNA Seq aims to represent a single condition, technology, or species and to discover cellular phenotypes. This enables the systematic reconstruction of cellular taxonomies in the human body. The biggest challenge is to identify subpopulations from multiple datasets. It is difficult to distinguish between changes in the composition of cell types in a sample and the expression changes within a given cell type when analyzing multiple datasets at the same time. Therefore, powerful new methods and a computational strategy are needed for learning between multiple datasets.
For example, zero-inflated differential expression tests have been tailored to scRNA-seq data to identify changes within a single cell type and clustering approaches detect proportional shifts across conditions if cell types are conserved. The methods should make it possible to learn between several data sets at the same time to facilitate a comparative analysis afterwards. Multivariate methods are used here, for example. One can identify gene correlation patterns that are conserved across data sets and embed cells in a common low-dimensional space (e.g., through CCA). CCA enables information from several data sets to be displayed consistently (linear combination of features in data sets that are maximally correlated). Data sets are treated as multiple measurements of a gene-gene covariance structure, and one looks for patterns that are common to the data sets. Multi-Set CCA enables the integration of several data sets.
scRNA Seq data are generally noisier and more complex than bulk RNA Seq data and therefore computational more challenging. 


### e.	What problem does the workflow at hand address (the Seurat vignette linked above)?

The Seurat workflow analyzes a dataset of Peripheral Blood Mononuclear Cells(PBMC). Raw single-cell expression data are used as an input. Aim of the workflow is finding clusters in the data with a graph-based clustering method. It uses a combination of feature selection, dimensionality reduction and clustering algorithms to identify cell types. 
We use the Seurat workflow to robustly separate different cell types in the sample. The target is to integrate the whole workflow from data pre-processing to separation of cell types and the assignment of cell types to the clusters. It used highly variable features to get reliable results that can be compared between different samples.

## 2.	Replication. To replicate the tutorial, you need to reproduce all figures presented in the workflow. Address at least the following questions

### a)	Is a replication of the tutorial possible? Compare the tutorial against the rules/recommendations from Sandve et al. 2013.; comment on the clarity of the description and documentation.

#### 10 rules/recommendations from Sandve et al. 2013:

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
The way the Seurat packages deals with randomness, was confusing for us. The seed is not set in the notebook, but hidden in the Seurat package. This leeds to the strange situation that for example the Jack Straw Plot differs visibly, if not significantly, from the tutorial, but there is no random change when run again as expected. The difference is due to the different handling of seed values on different operating systems.
For us it had been nicer to set the seed in the final notebook to avoid confusion.

7. Always Store Raw Data behind Plots
Most of the plots just require the Seurat object, which is covered in Rule 5. It would be nice to store additional data, like used features and cluster labels' in an easier to read form. In the current form they are in the code. So everything is reproducible, but not always easy to find.

8. Generate Hierarchical Analysis Output, Allowing Layers of Increasing Detail to Be Inspected
The Tutorial is more a way to show the capabilities of the package than tha good way to represent results. Therefore rule is not applied.

9. Connect Textual Statements to Underlying Results
Similar to rule 8 this is missing as a thorough explanation of the results is not part of the tutorial.

10. Provide Public Access to Scripts, Runs, and Results
The code is available in the GitHub repository.

### b)	How did you set up the required environment?
First we set up the environment. While there are different options for this we decided to use anaconda (miniconda) for this task. We set up the conda forge channel as default and run the following command:

`mamba create -n seurat r-base=4.1.0 r-ggplot2=3.3.5 patchwork=1.1.1 SeuratObject=4.0.2 Seurat=4.0.4 dplyr=1.0.7 r-knitr=1.33`

to create the environment. We kept the package versions used in the tutorial to avoid the risk of incompatibilities and tried to keep the environment as small as possible. The environment was exported as .yaml file to allow easy replication.

mamba is used instead of conda to to speed up the environment management. This is especially relevant as we use the large conda forge channel.

It is also required to create the folder structure manually. The notebooks have to be in a subfolder (scripts in our case), the exact location of the data files has to be set in the notebook and an output folder must be present before running the code.

### c)	Explain all the steps of the vignette in your own words.

### Configuration steps

First we set up the environment. While there are different options for this we decided to use anaconda (miniconda) for this task. We set up the conda forge channel as default and run the following command

`mamba create -n seurat r-base=4.1.0 r-ggplot2=3.3.5 patchwork=1.1.1 SeuratObject=4.0.2 Seurat=4.0.4 dplyr=1.0.7 r-knitr=1.33`

to create the environment. We kept the package versions used in the tutorial to avoid the risk of incompatibilities and tried to keep the environment as small as possible. The environment was exported as .yaml file to allow easy replication.

Mamba is used instead of conda to to speed up the environment management. This is especially relevant as we use the large conda forge channel.

It is also required to create the folder structure manually. The notebooks have to be in a subfolder (scripts in our case), the exact location of the data files has to be set in the notebook and an output folder must be present before running the code.

#### Work flow description (basic details in R Notebook)

The workflow runs in the following basic steps: 

##### Pre-processing
First we discard low quality cells, e.g. cells with very few genes, or suspected doublets and triplets with an high amount of genes. We normalize the data using a logarithmic method. 

#### Identification of highly variable features and PCA

Next we limit the further steps to the features showing a high variability. This ensures that we get less noise from the other features and makes the later calculation faster and more importantly leads to better results.

For the PCA the data from the variable features is normalized (mean = 0, variance = 1). After the PCA is performed the results are visualized. The loadings are presented, a scatter plot of the first PCs is created and heatmaps for the first PCs are created. This helps to identify where the relevant information can be found and which PCs and which loadings/genes are used later.

#### determine 'dimensionality' of dataset

After the PCA we need to know how many dimensions are needed for the further analyses. In the Tutorial 2 heuristic methods are used. The first is the Jack Straw plot and the second the elbow plot.
  * describe how the plots are used

#### Cluster the cells

A graph based clustering approach is used for the identified relevant PCs. As we do not know a priory how many clusters/cell type we will find this is a sensible approach.
The clusters are displayed using the UMAP non linear dimension reduction.

#### finding differentially expressed features

Finally we want to identify the clusters. We use VlnPlot to see in which cluster the genes with the highest loadings are expressed. Additionally the genes with the highest fold change for each cluster are calculated. With additional information about the cell types expected in the sample the clusters can be assigned to cell types.

#### Assess reproducibility of code

cite Sandve 2013
The rules form the Article are used to assess how easy the code can be reproduced and reused.

## Results

???? Plots oder Figures können wir hier nicht reinstellen oder??????

### das von meinem Word > nachdem du einiges geändert hast und ich da noch nicht ganz durchblicke würde ich vorschlagen du ergänzt da noch was du magst zu dem was du geschrieben hast...und den Rest geben wir weg

The Seurat workflow analyses a dataset of Peripheral Blood Mononuclear Cells (PBMCs) from 10X Genomics and consists out of 2700 single cell data sequenced with Illumina's NextSeq 500.
Setting up the Seurat objects begins with reading in the data with the data from the Cell ranger pipeline from 10 genomics. CellRanger is an analysis pipeline that processes single-cell data to align reads, generate feature barcode matrices, performs clustering and other secondary analysis. The output is an UMI (Unique Molecular Identifier) count matrix. UMIs are short sequences that are added to DNA sequences to uniquely identify the DNA molecule. The values in the matrix represent the number of molecules for each feature (gene, row) detected in each cell (column). The Count Matrix is used to create the Seurat object, which contains both the data and the analyses (e.g., PCA, clustering results), for a data set.
The preprocessing step is characterized by the selection and filtering of cells based on QC metrics, normalization of data and scaling or detection of highly variable features. QC metrics are defined by the user and are, for example: the number of unique genes in each cell (cells with poor quality or empty droplets often contain only a few genes, duplicates or multiplets can have very high gene counts). Unique genes correlate strongly with the number of molecules detected within a cell. Another metric is the percentage of reads that maps to the mitochondrial genome. Dying cells and cells with poor quality often have extensive mitochondrial contamination. The calculation is done with PercentageFeatureSet(), a function that calculates the percentage of counts from a given feature set. (All genes starting with MT).
Feature-feature relationships are visualized then by Feature Scatter and subsets are defined with specific filter criteria. After the filter steps, the data must be normalized. By default, the "LogNormalize" method is used for this, which normalizes the feature expression measurements for each cell by the total expression, multiplies by the scale factor (default 10000) and log-transforms the result.
A subset of features is defined with high cell-to-cell variation in the expression (highly expressed in some cells and lowly expressed in others). The further focus of the analysis is then placed on these genes. Therefore, the FindVariablesFeatures () function is used which allows directly modeling of the mean-variance relationship in single-cell data. Default is 2000 features per dataset, used e.g., for PCA. Scaling of the data is performed afterwards.
Scaling is the linear transformation which is the standard before dimension reduction techniques like PCA (mean expression across cells is 0 and variance across cells is 1). Scaling is intended to achieve equal weighting in downstream analysis so that highly expressed genes do not dominate.
The next step is to perform linear dimensions reduction (PCA) on the scaled data. By default, the previously defined features are used, but a new subset can also be defined. Cell and features that define PCA are visualized with DimHeatmap. DimHeatmap() allows visualization of the primary source of heterogeneity in the dataset and can help to decide which PCs to use for further downstream analysis. Cells and features are ordered by their PCA scores, which also increases speed, which is especially important for large datasets.
To eliminate technical noise in the individual features and to determine the dimensionality of the dataset, Seurat clusters the cells based on their PCA scores, each PC representing a so-called "metafeature" which combines information about a correlated feature set. In a resampling test, a subset of the data (1% by default) is randomly selected and the PCA " is repeated. A "null distribution of feature scores is identified...the procedure is repeated...this procedure allows to identify significant PCs that contain many low p-values.
Visualization allows the comparison of the distribution of p-values for each PC with a unit distribution. Significant PCs show a strong enrichment of features with low p-values (solid curve above the dashed line), a sharp drop-off in significance after the first 10 - 12 PCs can be seen. The elbow plot shows a ranking of PCs based on the percentage of variance explained (cutoff between PC 7 -12)
To identify the true dimensionality of a dataset is a challenge. Therefore, several approaches should always be included:
1. Supervised (exploring PCs to determine sources for heterogeneity), e.g., in combination with GSEA
2. Statistical test based on the random null model (time consuming for large datasets, no clear PC cutoff)
3. Heuristic (calculated instantly)
10 PCs were chosen for the cutoff, the following still needs to be considered: Dendritic cell and NK aficionados may recognize that genes strongly associated with PCs 12 and 13 define rare immune subsets (e.g., MZB1 is a marker for plasmacytoid DCs). However, these groups are so rare that they are difficult to distinguish from background noise for a dataset of this size without prior knowledge.
Downstream analysis is then repeated with a different number of PCs (10, 15, or even 50!). For clustering the cells Seurat uses a graph-based clustering approach. The distance metric remains the same and is based on the previously identified PCs. The approach is based on manuscripts applied to graph-based clustering approaches with scRNA seq data and is characterized by a graph structure, e.g., KNN nearest neighbors. (Edges drawn between cells with similar feature expression patterns) and then graph is partitioned into highly interconnected communities. 
First, a KNN graph is constructed based on the Euclidean distance (in PCA space), edge weights are refined between all two cells based on the overlap in their neighborhood. This step is done using the FindNeighbors() function with the first 10 PCs as input. 
To cluster the cells, modularity optimization techniques such as the Louvain algorithm or SLM are applied. Groups are iteratively clustered with the goal of optimizing the standard modularity function. The FindClusters() function is used for this and sets the granularity for downstream clustering (higher values, more clusters). Settings: 0.4 - 1.2 for single-cell datasets with about 3K cells. The larger the dataset the more likely is an optimal solution. With the Idents () function the clusters can be found.
Nonlinear dimensional reduction is performed with UMAP/tSNE. The goal of nonlinear dimensional reduction techniques is to learn underlying manifold of the data and place similar cells together in low-dimensional space. Cells in the previously determined clusters are to be displayed co-localized in the plots. The same PCs are used as input as for the cluster analysis.
Then, biomarkers must be found which define clusters through differential gene expression. By default, positive and negative markers of a single cluster are identified compared to all other cells. Therefore, a threshold value must be set: A trait must be recognized to a minimum percentage in one of the two cell groups and must be expressed differently between the two groups. (Both set to 0 leads to dramatic time savings). The maximum number of cells can also be defined which increases the speed.
Finally, visualization of the marker expression is performed: The violin plot shows the expression probability distribution over the clusters. The feature plot visualizes the feature expression in a tSNE or PCA plot. RidgePlot (), CellScatter () and DotPlot () can be further options. DoHeatmap () can be used to generate an expression heatmap for cells and features. The top 20 markers for each cluster are plotted then. 
As the last step in the workflow, canonical markers are sent to assign the unbiased clustering to known cell types (9 assignments: "Naive CD4 T", "CD14 + Mono", "Memory CD4 T", "B", "CD8 T", " FCGR3A + Mono ", "NK", "DC", "Platelet"). The new clusters are visualized as umap with DimPlot ().

## 3.	Expanding the work. Find a publicly available data set and apply the same workflow. You may need to adapt some of the code to make it work.

### Analysis on own dataset

A new dataset was introduced to test the possibility to reuse the code.

## Results
The dataset we used is from 10X Genomics. 
https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-2-standard-2-1-0

Description:
Cells from a combined cortex, hippocampus and sub ventricular zone of an E18 mouse.
•	931 cells detected
•	Sequenced on Illumina HiSeq2500 with approximately 56,000 reads per cell
•	26bp read1 (16bp Chromium barcode and 10bp UMI), 98bp read2 (transcript), and 8bp I7 sample barcode

### a)	What challenges did you face when applying the workflow to a new data set?

One challenge was certainly to define the features/markers for the final analysis. In the original tutorial, this point was not explicitly described or the markers were already predefined in the code. When using the new dataset, it was of course necessary to deal with it more intensively. For this purpose it was necessary to look at the code in more detail in order to be able to make a meaningful selection. 

### b)	What code modifications were required?
* First the path had to be changed to the folder with the new dataset
* In the new dataset the mitochondrial genes start with "mt-". The lower case is not recognized by the original pattern. The code was adapted to find the desired genes in the    new dataset
* The new dataset has a higher number of features. Therefore the filter was redefined according to the plots 
* The name of the saved Seurat object is changed to seperate the different datasets
* Assigning marker genes to the clusters: In this step we do not have enough background information to link our marker genes with cell types. Therfore we mark the clusters with the best marker gene.

### c)	Are the results comparable to the results of the original tutorial, or do they deviate in some
### unexpected ways?

The replication of the tutorial with the new dataset yielded results with the new dataset which unfortunately were not quite as nice as in the original workflow. Only PC 1 showed a good separation in the new dataset, which made the plots more difficult to interpret and the results not quite as impressive. 

### d)	Discuss all the results and interpret them

It is possible to reproduce the tutorial and use the method on a different dataset. While it worked well there are some options to make it easier. First the form of a notebook makes it difficult to see where parameters are set and manual interventions are necessary. It is hard to make changes without checking the whole code. In our opinion it would be easier to write the code in functions and use all variables as function parameters.

It would be possible to do the analysis in 2 well defined steps. Step 1 starting from pre-processing and ending with the analysis of the PCs. At this point manual intervention is done. The rest of the analyses can be performed without manual intervention in a second function or method.

### Literature

Luecken and Theis 2019. Current best practices in single-cell RNA-seq analysis: a tutorial. Molecular
Systems Biology. 15(6):e8746. doi:10.15252/msb.20188746.

Butler et al. 2018. Integrating single-cell transcriptomic data across different conditions, technologies,
and species. Nat Biotechnol. 36(5):411–420. doi:10.1038/nbt.4096.

Sandve et al. 2013. Ten Simple Rules for Reproducible Computational Research. PLOS Computational
Biology. 9(10):e1003285. doi:10.1371/journal.pcbi.1003285.
