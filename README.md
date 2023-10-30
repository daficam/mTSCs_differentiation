# mTSCs Differentiation timecourse ::: Dropseq

> Analysis by:: 'Malwina Prater'
> Contact E-mail:: 'mn367@cam.ac.uk'
> Application Type:: 'scRNA-seq'
> Project Owners:: 'Malwina Prater, Dafina Angelova, Russell Hamilton & Steve Charnock-Jones'
> Single Cell Platform:: 'Dropseq '
> Sequencing Platform:: 'Illumina '
> Sequencing Setup:: 'variable between batches'



# Methods


## Bioinformatics Methods: Processing of scRNA-Seq

Reads trimmed to 50bp length before aligning. This improved quality of data and increased number of genes per cell (mean nFeature_RNA = 810, median nFeature_RNA = 614 ).


## Bioinformatics Methods: processing

Samples with less than 40% alignment rate were removed.  Expression of all cells from each sample was added to create bulk-like RNA-seq dataset. Cumulative sample counts underwent variance stabilizing transformation (function `vst` from R package "DESeq2") and PCA was computed (prcomp function from R package "stats", v3.5.3). Sample to sample correlation was calculated by using rcorr function on transformed counts (R package hmisc v4.5.0) to compute a matrix of Pearson's r correlation coefficients for all possible pairs of columns of a matrix. Sample replicates were merged after confirmation of labbook notes using sample correlation and number of shared UMIs (r2 > 0.85 & number of shared UMIs > 100). scRNA-seq analysis was performed using package "Seurat" (v.3.2.0). We removed cells with fewer than 500 genes detected, leaving 136291 single cells in total. For all single-cell analysis, we performed the same initial normalization using R package “sctransform” (v.0.3.2)  that performs regularized negative binomial regression for the UMI counts (ref). Sctransform also calculated scaled expression for downstream dimensional reduction. Dimensional reduction (PCA, UMAP)  and finding neighbouring cells was performed with Seurat using the first 30 dimensions. To run UMAPs parameters these were used: reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 25, repulsion.strength = 2, spread = 2L. 

>Outlier detection:  

Mean expression of genes per sample was calculated for quality control and 1 sample outlier was detected on sc-to-bulk PCA (1 sample outlier), as well as single-cell PCA and UMAP. Additional outlier was identified during assessment of heterogeneity of T0 cell fraction, and was removed. 

>Batch integration:  

Since the dataset was sequenced in batches, batch effect was investigated and identified. Several tools were compared for best performance in batch correction/integration: Seurat Integration (with/without reference), Harmony, scaling out and BUSSeq. Harmony (v.1.0) was best to correct for batch effect with settings: reduction = "pca", assay.use="SCT". 

>Resolution:  

Resolution of clustering was assessed using R package clustree (v.0.4.3), separation on UMAP and number of differentially expressed genes in pairwise clusters to avoid under and overclustering. Resolution 0.8 was chosen for Together dataset as well as Inhibit only and Removal only, while resolution 0.4 was chosen for T24 and T0 only datasets. Clusters smaller than 100 cells were removed from analysis (cluster 20 in Together dataset). 

>Cell Cycle Scoring: 

Cells were scored for the cell cycle phase using the CellCycleScoring function in Seurat. R Package Tempora (v.0.1.0) was used to estimate cluster order. This cluster order was used in all subsequent figures that contain cells from different time points.  

>Pseudotime:  

Pseudotime trajectory was calculated using R package Monocle3 (v.0.1.2) using top 3000 most variable genes, with Msonocle object using UMAP embedding and clusters calculated in Seurat for consistency.   

>DEG and Term Enrichment analysis: 

Cluster markers were calculated for all cluster pairs using Seurat function “FindMarkers” with default log2 fold threshold 0.25. GO analysis and pathway enrichment was performed using enrichR R package (v.3.0) and databases: "KEGG_2019_Mouse", "GO_Biological_Process_2018", "GO_Cellular_Component_2018",  "GO_Molecular_Function_2018",  "BioCarta_2016",  "WikiPathways_2019_Mouse". Heatmaps were generated using R package ComplexHeatmap (v.1.20.0). Terms for reduced heatmaps were chosen on the basis of relevance to processes related to stem cell differentiation, mitosis, transcription, translation, cell remodelling and mitochondrial respiration. In GO databases enriched pathways were often related to the same process. They were grouped by parent term using rrvgo R package (threshold of reduction : 0.7) and reduced to child terms that were enriched in the most clusters. 

To understand cluster-specific related pathways, for each cluster and for all cells, average normalized expression of each gene was calculated to produce scatter plots comparing each cluster with average normalized expression of all dataset. Normalized average counts were log2 transformed. A cut-off >= 1 of abs(log2 mean expression Individual Cluster  - log2 mean total expression)  was used to identify variable genes, and additional filtering was applied to include only genes that were identified as differentially expressed in any pairwise comparison between clusters. For each cluster, genes with elevated and lower expression levels compared to the mean expression of the total dataset were separately used for pathway enrichment (WikiPathways_2019_Mouse). For rare events when both elevated and lower expressed genes were associated with the same pathway, either elevated or lower expressed group was chosen on the basis of lower FDR. 

Transcription Factor regulons and gene regulatory networks were inferred using R package SCENIC (v.1.1.2-2) and GRNBoost2 algorithm from python library Arboreto (v.0.1.6). This method (GRNBoost2 algorithm) was chosen on the basis of performance, time and memory usage and low sensitivity to the presence of dropouts (Pratapa et al, 2020. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. Nature Methods Vol 17, p147-154).
