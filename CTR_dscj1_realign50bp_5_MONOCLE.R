#!/storage/Software/packages/anaconda3/bin/Rscript

message("--------------------------------------------------------------------------------")
message("+                         MONOCLE on SCT dataset                                ")
message("+-------------------------------------------------------------------------------")

rm(list=ls())
gc()
options(bitmapType='cairo')


library(Seurat)
library(ComplexHeatmap)
library(circlize)
library("dplyr")
library("methods")
library("utils")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Matrix")
library("matrixStats")
library("Seurat")  
library("useful")
library("reshape2")
library("DESeq2")   
library("ggalt")
library("Matrix.utils")
library("scran")
library("scater")
library("monocle3")


baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
Project <- "CTR_dscj1_MONOCLE"
resDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/MONOCLE_SCT"
setwd(resDir)

col_expt <- c("0.0"="black", "I.1"="yellow", "R.1"="yellow2", "I.4"="orangered", "R.4"="red", "I.24"="lightblue2", "R.24"="lightblue4", "I.36"="blue", "R.36"="darkblue","I.48"="purple", "R.48"="orchid")  


#dataset <- "harmony.orig_500g"
#dataset <- "harmony_R_500g"
dataset <- "harmony_I_500g"




message("--------------------------------------------------------------------------------")
message("+                                Load objects:::                                ")
message("+-------------------------------------------------------------------------------")
reso <- "SCT_snn_res.0.8"

#matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")
#matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_orig/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729.Rds")

if (dataset == "harmony.orig_500g"){
  dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")
  Idents(dropseq.integrated) <- dropseq.integrated@meta.data[[reso]]
} else if (dataset == "harmony_R_500g"){
  dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_R.Rds")
  Idents(dropseq.integrated) <- dropseq.integrated@meta.data[[reso]]
  # from tempora for R 500g::: 
  Idents(dropseq.integrated) <- factor(Idents(dropseq.integrated), levels= c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8"))
} else if (dataset == "harmony_I_500g"){
  dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_I.Rds")
  Idents(dropseq.integrated) <- dropseq.integrated@meta.data[[reso]]
  # from tempora for I 500g::: 
  Idents(dropseq.integrated) <- factor(Idents(dropseq.integrated), levels= c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0"))
} else { "You didn't provide valid dataset!!"}





### Assign UMAP coordinate
UMAP_embeddings <- dropseq.integrated@reductions[["umap"]]@cell.embeddings
head(UMAP_embeddings) # dim:  123443      2
saveRDS(UMAP_embeddings, paste(Project, dataset, "UMAP_embeddings.Rds", sep = "_"))


sampleTable <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/sampleTable_MERGED_53.csv")
metadata <- dropseq.integrated@meta.data
saveRDS(metadata, paste0(Project, dataset, "_metadata.Rds"))

variable_genes_3000 <- VariableFeatures(dropseq.integrated)
gene_annotation <- as.data.frame(variable_genes_3000, row.names = variable_genes_3000)
colnames(gene_annotation) <- "gene_short_name"
dim(gene_annotation) # 5000 genes annotated as in PCA

expression_matrix <- GetAssayData(dropseq.integrated, assay = "SCT", slot = "data")
expression_matrix <- expression_matrix[rownames(gene_annotation), ] 



cell_metadata <- as.data.frame(dropseq.integrated@assays[["SCT"]]@counts@Dimnames[[2]], row.names = dropseq.integrated@assays[["SCT"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
dim(cell_metadata) # [1]   127781      1/6
head(cell_metadata)


#saveRDS(gene_annotation, paste0(Project, dataset, "_gene_annotation_3k.Rds"))
#saveRDS(cell_metadata, paste0(Project, dataset, "_cell_metadata.Rds"))
#saveRDS(expression_matrix, paste0(Project, dataset, "_expression_matrix_3k.Rds"))


#expression_matrix_all <- readRDS(paste0(Project, dataset, "_expression_matrix.Rds"))
expression_matrix <- readRDS(paste0( Project, dataset, "_expression_matrix_3k.Rds"))
cell_metadata <- readRDS(paste0(Project, dataset, "_cell_metadata.Rds"))
#expression_matrix <- readRDS(paste0(Project, dataset, "_expression_matrix.Rds"))
gene_annotation <- readRDS(paste0(Project, dataset, "_gene_annotation_3k.Rds"))
all_metadata <-readRDS(paste0(Project, dataset, "_metadata.Rds"))

#filtered.counts <- expression_matrix_all[ Matrix::rowSums(expression_matrix_all>0)>50,]
#dim(filtered.counts)

#gene_annotation_filt <- as.data.frame(rownames(filtered.counts))
#colnames(gene_annotation_filt) <- "gene_short_name"
#rownames(gene_annotation_filt) <- gene_annotation_filt$gene_short_name

# remove cells from cluster 20!
head(all_metadata,2)
all_metadata_f <- all_metadata[all_metadata$SCT_snn_res.0.8 != "20",]
cell_metadata_f <- subset(cell_metadata, cell_metadata$barcode %in% rownames(all_metadata_f))
dim(cell_metadata)
filtered.counts_f <- filtered.counts[, colnames(filtered.counts) %in% rownames(all_metadata_f)]
dim(filtered.counts_f)

#cds_from_seurat_f     <- monocle3::new_cell_data_set(filtered.counts_f, 
#                                                   cell_metadata = cell_metadata_f,
#                                                   gene_metadata = gene_annotation_filt)

cds_from_seurat     <- monocle3::new_cell_data_set(expression_matrix, 
                                                     cell_metadata = cell_metadata,
                                                     gene_metadata = gene_annotation)


head(cds_from_seurat@colData)

cds_from_seurat@colData$Sample      <- rownames(cds_from_seurat@colData)
cds_from_seurat@colData$Sample      <- gsub("_[ACGTN]{12}", "", cds_from_seurat@colData$Sample)
cds_from_seurat@colData$Batch       <- sampleTable[match(cds_from_seurat@colData$Sample, sampleTable$sampleLabels2), ]$Merged_batch
cds_from_seurat@colData$Experiment  <- sampleTable[match(cds_from_seurat@colData$Sample, sampleTable$sampleLabels2), ]$CellType
cds_from_seurat@colData$Phase       <- all_metadata[match( rownames(cds_from_seurat@colData), rownames(all_metadata)), ]$Phase
cds_from_seurat@colData$SCT_snn_res.0.8 <- all_metadata[match( rownames(cds_from_seurat@colData), rownames(all_metadata)), ]$SCT_snn_res.0.8
#saveRDS(cds_from_seurat, paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_filt15k_noC20.Rds"))
#saveRDS(cds_from_seurat, paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_3k_noC20.Rds"))
#saveRDS(cds_from_seurat, paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_3k_noC20.Rds"))
#saveRDS(cds_from_seurat, paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_3k.Rds"))

rm(matrix.su)




#cds_from_seurat <- readRDS(paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_3k_noC20.Rds"))
#cds_from_seurat <- readRDS(paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_filt15k.Rds"))
#cds_from_seurat <- readRDS(paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_filt15k_noC20.Rds"))
# for R or I:
#cds_from_seurat <- readRDS(paste0(Project, "_", dataset, "_", "cds_from_seurat_ini_3k.Rds"))




set.seed(42)
gc()

## Step 1: Normalize and pre-process the data
cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 50, scaling = TRUE)
#saveRDS(cds_from_seurat, paste0(Project, "_", dataset, "_", "cds_from_seurat_preprocessed_3k.Rds"))
gc()

## Step 2: Reduce the dimensions using UMAP
cds_from_seurat <- reduce_dimension(cds_from_seurat, reduction_method = "UMAP", preprocess_method = "PCA", umap.min_dist = 0.01, verbose = TRUE, umap.n_neighbors = 30L)
#saveRDS(cds_from_seurat, file ="Monocle3_cds_from_seurat_redDim_3k.rds")
gc()

## Step 3: Cluster the cells
cds_from_seurat <- cluster_cells(cds_from_seurat, k = 100, verbose = TRUE) # num_iter = 100, 
#saveRDS(cds_from_seurat, file ="Monocle3_cds_from_seurat_UMAP_clust_k100_3k.rds")




#### issue with partitions... there is a messy way around:::

unique(cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]]) 
table(cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]]) 


# you can plot to check which cells are in whihc partition
pData(cds_from_seurat)$Partition <- cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]]
plt_partition_cds_from_seurat <- plot_cells(cds_from_seurat, color_cells_by = "Partition", group_cells_by = "partition", label_leaves=FALSE, label_roots = FALSE, label_cell_groups = TRUE, label_branch_points=FALSE, reduction_method = "UMAP", group_label_size = 6, alpha = 0.3)

pdf(paste("Monocle3_UMAP__cds_partition", dataset, "3k", ".pdf", sep="_"), width=15,height=10)
par(bg=NA)
plt_partition_cds_from_seurat
dev.off()


pData(cds_from_seurat)$Experiment <- factor(pData(cds_from_seurat)$Experiment, levels = c( "0.0", "R.1", "R.4", "R.24", "R.36", "R.48","I.1", "I.4", "I.24", "I.36", "I.48"))
plt_expt_cds_from_seurat <- plot_cells(cds_from_seurat , color_cells_by = "Experiment", group_cells_by = "cluster", label_leaves=FALSE, label_roots = FALSE, label_cell_groups = TRUE, label_branch_points=FALSE, reduction_method = "UMAP", group_label_size = 6, alpha = 0.3)

plt_expt_cds_from_seurat     <- plt_expt_cds_from_seurat + scale_colour_manual(values = col_expt)


pdf(paste("Monocle3_UMAP__cds_Expt", dataset, "3k", ".pdf", sep="_"), width=15,height=10)
par(bg=NA)
plt_expt_cds_from_seurat
dev.off()





if (dataset == "harmony.orig_500g"){
  # for harmony orig 500g reso 0.8 (all clusters, nothing removed):::
  Col_cluster <- c( "#FFFF99" , "#d9f0a3" , "#bfd3e6" , "#fa9fb5" ,  "#88419d", 
                    "#0868ac" ,  "#f768a1" ,  "#4d004b" ,  "#00B3FFFF", "#FF00E6FF",
                    "#FFE500FF", "#004529"  , "#78c679" ,  "#FF0000FF", "#810f7c" ,
                    "#081d58" ,  "#00FFB2FF", "#fde0dd" , "#FF9900FF" ,   "#7F00FFFF"  )
  
  pData(cds_from_seurat)$Clust <- factor(pData(cds_from_seurat)$SCT_snn_res.0.8, levels = c( "0", "1", "2", "3", "4", "5","6", "7", "8", "9", "10", "11", "12","13","14", "15", "16", "17", "18", "19"))
} else if (dataset == "harmony_R_500g"){
  Col_cluster<-c( "#8c96c6" , "#00B3FFFF", "#4d004b" ,  "#238443" , "#f768a1"  , 
                  "#88419d"  , "orange"  , "#fde0dd"  , "yellow2" , "#d4b9da",
                  "#78c679"  , "#d9f0a3" ,  "#4eb3d3" , "#081d58" ,  "#00FFB2FF",
                  "#0868ac"  )
  pData(cds_from_seurat)$Clust <- factor(pData(cds_from_seurat)$SCT_snn_res.0.8, levels = c( "0", "1", "2", "3", "4", "5","6", "7", "8", "9", "10", "11", "12","13","14", "15"))
  
} else if (dataset == "harmony_I_500g"){
  # for res 0.8 I only ::: 500g
  Col_cluster<-c( "#8c96c6" , "#00B3FFFF", "#4d004b" ,  "#238443" , "#f768a1"  , 
                  "#88419d"  , "orange"  , "#fde0dd"  , "yellow2" , "#d4b9da",
                  "#78c679"  , "#d9f0a3" ,  "#4eb3d3" , "#081d58" ,  "#00FFB2FF",
                  "#0868ac" ,  "#fa9fb5" ) # ,  "#FF0000FF", "#00FF19FF"
  pData(cds_from_seurat)$Clust <- factor(pData(cds_from_seurat)$SCT_snn_res.0.8, levels = c( "0", "1", "2", "3", "4", "5","6", "7", "8", "9", "10", "11", "12","13","14", "15", "16"))
  
  
} else { "You didn't provide valid dataset!!"}

names(Col_cluster) <- levels(pData(cds_from_seurat)$Clust)





plt_Clust_cds_from_seurat <- plot_cells(cds_from_seurat , color_cells_by = "Clust", group_cells_by = "cluster", label_leaves=FALSE, label_roots = FALSE, label_cell_groups = TRUE, label_branch_points=FALSE, reduction_method = "UMAP", group_label_size = 6, alpha = 0.3)

plt_Clust_cds_from_seurat     <- plt_Clust_cds_from_seurat + scale_colour_manual(values = Col_cluster)


pdf(paste("Monocle3_UMAP__cds_SeuratClust", dataset, "3k", ".pdf", sep="_"), width=15,height=10)
par(bg=NA)
plt_Clust_cds_from_seurat
dev.off()





levels(cds_from_seurat@clusters$UMAP$partitions)
cds_from_seurat@clusters$UMAP$partitions[cds_from_seurat@clusters$UMAP$partitions == "2"] <- "1"
cds_from_seurat@clusters$UMAP$partitions[cds_from_seurat@clusters$UMAP$partitions == "3"] <- "1"
unique(cds_from_seurat@clusters$UMAP$partitions)
cds_from_seurat@clusters$UMAP$partitions <- factor(cds_from_seurat@clusters$UMAP$partitions, levels = c("1"))


message("--------------------------------------------------------------------------------")
message("+                         Learn graph                                           ")
message("+-------------------------------------------------------------------------------")

cds_from_seurat     <- learn_graph(cds_from_seurat,  use_partition = FALSE, close_loop = FALSE)
#saveRDS(cds_from_seurat, file ="Monocle3_cds_from_seurat_LG_UsePartF_LoopT__3k_noC20_rmPartitions.rds")
saveRDS(cds_from_seurat, file =paste("Monocle3_", dataset, "cds_from_seurat", "LG_UsePartF_LoopT__3k_rmPartitions.rds", sep = "_"))




# cds with monocle embedding"
#cds_from_seurat <- readRDS( "Monocle3_harmony.orig_500g_3kVarGenes/Monocle3_cds_from_seurat_LG_UsePartF_LoopT__3k_noC20_rmPartitions.rds")
#cds_from_seurat <- readRDS( "Monocle3_harmony_Inhibit_500g_3kVarGenes/Monocle3__harmony_I_500g_cds_from_seurat_LG_UsePartF_LoopT__3k_rmPartitions.rds")
#cds_from_seurat <- readRDS( "Monocle3_harmony_Remove_500g_3kVarGenes/Monocle3__harmony_R_500g_cds_from_seurat_LG_UsePartF_LoopT__3k_rmPartitions.rds")

# cds with SEURAT embedding"
#cds_from_seurat <- readRDS( "Monocle3_harmony.orig_500g_3kVarGenes/Monocle3_cds_mix___SeuratUMAP___LG_UsePartF_LoopT__3k_noC20.Rds")
#cds_from_seurat <- readRDS( "Monocle3_cds_mix___SeuratUMAP___LG_UsePartF_LoopFALSE__3k_noC20.Rds")
cds_from_seurat <- readRDS( "Monocle3_harmony_Inhibit_500g_3kVarGenes/Monocle3__I__cds_mix___SeuratUMAP___LG_UsePartF_LoopT__3k.Rds")
#cds_from_seurat <- readRDS( "Monocle3_harmony_Remove_500g_3kVarGenes/Monocle3__R__cds_mix___SeuratUMAP___LG_UsePartF_LoopT__3k_noC20.Rds")


pData(cds_from_seurat)$Clust <- pData(cds_from_seurat)$SCT_snn_res.0.8



#cds_from_seurat     <- learn_graph(cds_from_seurat,  use_partition = FALSE, close_loop = FALSE)
#saveRDS(cds_from_seurat, file ="Monocle3_cds_mix___SeuratUMAP___LG_UsePartF_LoopFALSE__3k_noC20.Rds")


message("--------------------------------------------------------------------------------")
message("+            calculate pseudotime - order cells             ")
message("+-------------------------------------------------------------------------------")

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="0.0"){
  cell_ids <- which(colData(cds)[, "Experiment"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds_from_seurat <- order_cells(cds_from_seurat, root_pr_nodes=get_earliest_principal_node(cds_from_seurat))


#saveRDS(cds_from_seurat, file =paste("Monocle3_", dataset, "cds_from_seurat", "LG_UsePartF_LoopF__3k_rmPart_noC20.rds", sep = "_"))



print("Plotting pseudotime")

plt_pseudotime     <- plot_cells(cds_from_seurat,     color_cells_by = "pseudotime",  label_leaves=FALSE, label_roots = FALSE, label_cell_groups = TRUE, label_branch_points=FALSE)


pdf(paste("Monocle3_UMAP__cds_PseudotimeTrajectory", dataset, "LoopFALSE", ".pdf", sep="_"), width=15,height=10)
par(bg=NA)
plot_cells(cds_from_seurat,  color_cells_by = "pseudotime",  label_leaves=FALSE, label_roots = FALSE, label_cell_groups = TRUE, label_branch_points=FALSE)
dev.off()

saveRDS(plt_pseudotime, paste("Monocle3_UMAP__cds_PseudotimeTrajectory", dataset, "x", ".Rds", sep="_"))


head(pData(cds_from_seurat))


Pseudotime <- pseudotime(cds_from_seurat, reduction_method = "UMAP")
saveRDS(Pseudotime, paste(Project, dataset, "Pseudotime_3k.Rds", sep = "_"))



message("--------------------------------------------------------------------------------")
message("+               calculate pseudotime - for clusters                             ")
message("+-------------------------------------------------------------------------------")


pData(cds_from_seurat)$Pseudotime        <- cds_from_seurat@principal_graph_aux[["UMAP"]]$pseudotime
head(pData(cds_from_seurat))
unique(pData(cds_from_seurat)$Clust)

clust_pstime_cds <- as.data.frame(pData(cds_from_seurat)[,c("Pseudotime", "Clust")])
clust_pstime_cds_agg <- aggregate(. ~ Clust, clust_pstime_cds, mean)
clust_pstime_cds_agg <- clust_pstime_cds_agg[order(clust_pstime_cds_agg$Pseudotime),]
cluster_order <- clust_pstime_cds_agg$Clust
as.character(cluster_order)

# Monocle cluster order on the basis of SEURAT embeddings for"::::: 
# R dataset:
# cluster_order <- c("14", "13", "2" , "3" , "7",  "0",  "15", "6" , "12", "11", "5" , "10", "9" , "4"  ,"8"  ,"1" ) 
# I dataset:
# cluster_order <- c("9",  "11" ,"16", "5" , "1" , "15", "6" , "12", "2" , "8" , "0" , "7" , "10", "14", "13", "4" , "3"  )
# Together dataset:
# cluster_order <- c( "12", "0" , "18", "15" ,"2" , "3" , "11", "17", "14" ,"7",  "9" , "10", "13" ,"8" , "1" , "19", "4",  "5" , "6" , "16" )
# together no loop:
#cluster_order <- c( "12" ,"0" , "18", "15", "3" , "2" , "11", "17", "9" , "13" ,"8" , "19", "5" , "7" , "4" , "6" , "16", "10", "14", "1") 


# Monocle cluster order on the basis of MONOCLE embeddings for"::::: 
# cluster order for R dataset:
#cluster_order <- c("14", "13", "2" , "15" ,"7",  "6" , "3" , "0" , "5"  ,"12", "10" ,"11", "9" , "4",  "1" , "8") 
# I dataset:
#  "9"  "11" "16" "5"  "1"  "15" "10" "12" "2"  "6"  "8"  "3"  "13" "14" "0"  "7"  "4"
# Together dataset:
# "12" "15" "0"  "18" "11" "17" "3"  "2"  "13" "7"  "14" "8"  "10" "1"  "9"  "5"  "16" "4"  "6"  "19"
saveRDS(cluster_order, paste(Project, dataset, "Pseudotime__cluster_order.Rds", sep = "_"))
write.csv(clust_pstime_cds_agg, paste(Project, dataset, "Pseudotime__clust_pstime_cds_agg.csv", sep = "_"))




message("--------------------------------------------------------------------------------")
message("+              Seurat embeddings for Monocle                                    ")
message("+-------------------------------------------------------------------------------")

# orig:
cds_from_seurat <- readRDS("Monocle3_harmony.orig_500g_3kVarGenes/Monocle3_cds_from_seurat_UMAP_clust_k100_3k_noC20.rds")
Seurat_UMAP <- readRDS("Monocle3_harmony.orig_500g_3kVarGenes/CTR_dscj1_MONOCLE_harmony.orig_500g_UMAP_embeddings.Rds")
Seurat_UMAP$UMAP_1 <- Seurat_UMAP$UMAP_1 * -1

# R:
cds_from_seurat <- readRDS( "Monocle3_harmony_Remove_500g_3kVarGenes/Monocle3__harmony_R_500g_cds_from_seurat_LG_UsePartF_LoopT__3k_rmPartitions.rds")
Seurat_UMAP <- readRDS("CTR_dscj1_MONOCLE_harmony_R_500g_UMAP_embeddings.Rds")

# I:
cds_from_seurat <- readRDS( "Monocle3_harmony_Inhibit_500g_3kVarGenes/Monocle3__harmony_I_500g_cds_from_seurat_LG_UsePartF_LoopT__3k_rmPartitions.rds")
Seurat_UMAP <- readRDS("CTR_dscj1_MONOCLE_harmony_I_500g_UMAP_embeddings.Rds")


# remove partitions before LG:
levels(cds_from_seurat@clusters$UMAP$partitions)
cds_from_seurat@clusters$UMAP$partitions[cds_from_seurat@clusters$UMAP$partitions == "2"] <- "1"
cds_from_seurat@clusters$UMAP$partitions[cds_from_seurat@clusters$UMAP$partitions == "3"] <- "1"
unique(cds_from_seurat@clusters$UMAP$partitions)
cds_from_seurat@clusters$UMAP$partitions <- factor(cds_from_seurat@clusters$UMAP$partitions, levels = c("1"))


cds_mix <- cds_from_seurat
Seurat_UMAP <- Seurat_UMAP[rownames(Seurat_UMAP) %in% rownames(pData(cds_mix)),]


head(cds_from_seurat@reducedDims@listData[["UMAP"]]) # UMAP from Seurat
head(cds_mix@reducedDims@listData[["UMAP"]]) # UMAP from Seurat
cds_mix@reducedDims@listData[["UMAP"]] <- Seurat_UMAP

### assign seurat clusters
list_cluster <- pData(cds_mix)$SCT_snn_res.0.8
names(list_cluster) <- rownames(pData(cds_mix))
head(list_cluster)
cds_mix@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
head(cds_mix@clusters@listData[["UMAP"]][["clusters"]])
head(cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]])

cds_mix  <- learn_graph(cds_mix,  use_partition = FALSE, close_loop = TRUE)
unique(cds_mix@clusters$UMAP$partitions)

# Together (orig):
#saveRDS(cds_mix, "Monocle3_harmony.orig_500g_3kVarGenes/Monocle3_cds_mix___SeuratUMAP___LG_UsePartF_LoopT__3k_noC20.Rds")

# R:
#saveRDS(cds_mix, "Monocle3_harmony_Remove_500g_3kVarGenes/Monocle3__R__cds_mix___SeuratUMAP___LG_UsePartF_LoopT__3k_noC20.Rds")

# I:
#saveRDS(cds_mix, "Monocle3_harmony_Inhibit_500g_3kVarGenes/Monocle3__I__cds_mix___SeuratUMAP___LG_UsePartF_LoopT__3k.Rds")










message("--------------------------------------------------------------------------------")
message("+           find genes that change in pseudotime             ")
message("+-------------------------------------------------------------------------------")

# read pre-calculated results:::
cds_pr_test_resSig <- readRDS(paste(Project, dataset,"prGraph_test_resSig.Rds", sep = "_"))
cds_pr_deg_ids <- readRDS( paste(Project, dataset, "cds_pr_deg_ids.Rds", sep = "_"))
gene_module_df <- readRDS( paste(Project, dataset,"gene_module.Rds", sep = "_"))




cds_pr_test_res = graph_test(cds_from_seurat, neighbor_graph="principal_graph", cores=4)
saveRDS( cds_pr_test_res,  paste(Project, dataset, "cds_pr_test_res.Rds", sep = "_"))

cds_pr_test_resSig <- subset(cds_pr_test_res, cds_pr_test_res$q_value < 0.05)
head(cds_pr_test_resSig)
saveRDS(cds_pr_test_resSig, paste(Project, dataset,"prGraph_test_resSig.Rds", sep = "_"))
write.csv(cds_pr_test_resSig, paste(Project, dataset,"prGraph_test_resSig.csv", sep = "_"))
cds_pr_test_resSig <- cds_pr_test_resSig[order(cds_pr_test_resSig$morans_test_statistic, decreasing = TRUE),]

cds_pr_deg_ids = row.names(subset(cds_pr_test_res, q_value < 0.05))
saveRDS(cds_pr_deg_ids, paste(Project, dataset, "cds_pr_deg_ids.Rds", sep = "_"))

plot_cells(cds_from_seurat, genes=c("Saa3", "Cdkn1c", "H2-D1", "Gjb2"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)



#  collect the trajectory-variable genes into modules:
#library(reticulate)
#source_python('louvain.py')

gene_module_df = monocle3:::find_gene_modules(cds_from_seurat[cds_pr_deg_ids,]) #, resolution=c(0,10^seq(-6,-1)), louvain_iter = 1)
write.csv(gene_module_df, paste(Project, dataset,"gene_module.csv", sep = "_"))
saveRDS(gene_module_df,  paste(Project, dataset,"gene_module.Rds", sep = "_"))


cell_group_df = tibble::tibble(cell=row.names(colData(cds_from_seurat)), cell_group=clusters(cds_from_seurat)[colnames(cds_from_seurat)])
rownames(cell_group_df) <- cell_group_df$cell
agg_mat = aggregate_gene_expression(cds_from_seurat, gene_module_df, cell_group_df)
module_dendro = hclust(dist(agg_mat))
#gene_module_df$module <- factor(module_dendro$module, levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_from_seurat,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


cell_group_df <- tibble::tibble(cell=row.names(colData(cds_from_seurat)), 
                                cell_group=colData(cds_from_seurat)$Clust)
agg_mat <- aggregate_gene_expression(cds_from_seurat, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")












head(colData(cds_from_seurat))
colData(cds_from_seurat)$Clust <- factor(colData(cds_from_seurat)$Clust, levels = cluster_order)



# check any interesting genes in module16:
gene_module_df <- gene_module_df[match(rownames(cds_pr_test_resSig) ,gene_module_df$id ),]
module13 <- gene_module_df[gene_module_df$module == 13,]$id
module16 <- gene_module_df[gene_module_df$module == 16,]$id
module2 <- gene_module_df[gene_module_df$module == 2,]$id
module5 <- gene_module_df[gene_module_df$module == 5,]$id
module3 <- gene_module_df[gene_module_df$module == 3,]$id
module15 <- gene_module_df[gene_module_df$module == 15,]$id
module12 <- gene_module_df[gene_module_df$module == 12,]$id

plot_cells(cds_from_seurat, genes=c("Ifit2", "Ifit1","Ifit3", 'Ifit3b' ), label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)

plot_cells(cds_from_seurat, genes=module2[1:12], label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)

plot_cells(cds_from_seurat, genes=module9[13:24], label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)

plot_cells(cds_from_seurat, genes=module9[25:36], label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)

plot_cells(cds_from_seurat, genes=module9[37:48], label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)



# some genes from module16 added to EPC progenitor lineage
EPC_genes <- c("Ascl2", "Gjb3","Tfap2c","Phlda2", "Tbfg1","Hbegf", "Tfrc", "Las1l", "Elf3", "Fgfbp1", "Sin3b", "Plac1", "Hmgn5", "Batf3", "Mbnl3", "Vav3","Ugp2","Prss8","Cxadr","Cmkir1","Ankrd2","Nacad","Pard6g","Serpinb9e","Nrk", "1600025M17Rik",        "Psap","Anxa4", "Slc22a18" ,"Slc16a1" , "Ctsb")


"Psap","Anxa4", "Slc22a18" ,"Slc16a1" , "Ctsb", "Lima1","Igsf8", "Reep6" ,"Stx3", 
#"Dnase1l3", "Ugp2", "Vgll3"  -Ugp2 - glycogen storage metabolism


plot_cells(cds_from_seurat, genes=EPC_genes, label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)

plot_cells(cds_from_seurat, genes=c("E2f8","Slc2a1", "Gjb2"), label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)



cds_sub <- choose_graph_segments(cds_from_seurat)


# For Inhibit:::
EPC_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% c("Las1l", "Hmgn5", "Hbegf","Plac1"),colData(cds_from_seurat)$Clust %in% c(9,16,11,1,15,6,8,0,3)]
# For Remove:::
EPC_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% c("Las1l", "Hmgn5", "Hbegf","Plac1"),colData(cds_from_seurat)$Clust %in% c(14, 2,0,6, 15,5,11,9,8,1,4,10)]


plot_genes_in_pseudotime(EPC_lineage_cds,
                         color_cells_by="Clust",
                         min_expr=0.2) + scale_colour_manual(values = Col_cluster)





LB_genes <- c("Gata2","Gcm1","Tmem37", "Bsg","Gjb2","Car4","Maged1","Wfdc2","F2rll1","P4hb", "Esx1","Tcfeb","Tead1","Tead4", "Tead5","Itm2a"  ,  "Tmem150a" ,"F2rl1","Apela","Rdh10" , "Tspan4","Brwd3", "Tbx3" ,"Uba6" , "Cdk5","Slc4a2","Kif5c","Fam129b","Gadd45g",    "Igfbp2", "Fhod3", "Col13a1")# "Dlx3","Lgals1","Fabp3","Ldha","Mtfr1l",  "Tspan4"  "Tbx3"    "Brwd3"

plot_cells(cds_from_seurat, genes=LB_genes,
           label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)


# For Inhibit:::
LB_lineage_cds <- choose_graph_segments(cds_from_seurat)
LB_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% c("Tmem37", "Gjb2", "Car4","F2rl1"),colData(cds_from_seurat)$Clust %in% c(9,11,16,1,15,6,12,2,7,4,13)]
# For Remove:::
LB_lineage_cds <- cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% c("Tmem37", "Gjb2", "Car4","F2rl1"),colData(cds_from_seurat)$Clust %in% c(14, 2,6, 15,11,9)]



plot_genes_in_pseudotime(LB_lineage_cds, 
                         color_cells_by="Clust", 
                         min_expr=0.2) + scale_colour_manual(values = Col_cluster)








plot_cells(cds_from_seurat, genes=c("Tmem37", "Gjb2", "Car4","F2rl1", "Las1l", "Hmgn5", "Hbegf","Plac1" ), label_cell_groups=FALSE,show_trajectory_graph=TRUE, label_leaves=FALSE, label_roots = FALSE, label_branch_points=FALSE)

















