#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()
options(bitmapType='cairo')


library(cowplot)
library(Seurat)
library(ggrepel)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)

message("--------------------------------------------------------------------------------")
message("+                  Load in dataset:::                                           ")
message("+-------------------------------------------------------------------------------")

dataset_name <- "Initial-no-merging"
matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")

dataset_name <- "Seurat_SCT_integration"
matrix.su <- readRDS("SCT_integration/dropseq.integrated_clust_CC.Rds")

dataset_name <- "Seurat_with_ref"
matrix.su <- readRDS("SCT_ref_0B0/dropseq.integrated_SCT_Ref0B0_clust.Rds")

dataset_name <- "Harmony"
matrix.su <- readRDS("HARMONY/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_clust.Rds")

dataset_name <- "Harmony_withRef"
matrix.su <- readRDS("HARMONY/harmony_theta4_ref0B0_matrix.su.Rds")

dataset_name <- "NonIntegr_ScaledOut"
matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformBatchScaledOut_clust.Rds")

dataset_name <- "BUSseq"
matrix.su <- readRDS("matrix.su_BUSseq_n.celltypes20.Rds")


message("--------------------------------------------------------------------------------")
message("+                  tables batches per clust                                     ")
message("+-------------------------------------------------------------------------------")

colnames(matrix.su@meta.data)

# if this not done, then do:::
#matrix.su <- FindNeighbors(matrix.su,  dims = 1:20) 
#matrix.su <- FindClusters(matrix.su, resolution = 1.25)
colnames(matrix.su@meta.data)
head(matrix.su@meta.data,3)


matrix.su@meta.data$CLUSTER <- matrix.su@meta.data$SCT_snn_res.1.25
matrix.su@meta.data$CLUSTER <- matrix.su@meta.data$integrated_snn_res.1.25
matrix.su@meta.data$CLUSTER <- matrix.su@meta.data$RNA_snn_res.1.25

tmp_df <- matrix.su@meta.data[,c("orig.ident" ,"Experiment", "Batch",  "CLUSTER")]

tmp_df$Experiment <- gsub( "\\.1", "\\.01" , tmp_df$Experiment)
tmp_df$Experiment <- gsub( "\\.4", "\\.04" , tmp_df$Experiment)
tmp_df$Experiment <- gsub( "\\.048", "\\.48" , tmp_df$Experiment)
tmp_df$Expt_batch <- paste(tmp_df$Experiment, tmp_df$Batch, sep = "_")
head(tmp_df)
tmp_df$CLUSTER <- as.numeric(as.character(tmp_df$CLUSTER))

tmp_df2 <- as.data.frame(table(tmp_df[,c("CLUSTER", "Batch")]))
tmp_df2 <- subset(tmp_df2, tmp_df2$Freq > 30)
tmp_df2 <- as.data.frame(table(tmp_df2$CLUSTER))
mean(tmp_df2$Freq)


message("--------------------------------------------------------------------------------")
message("+           Dispersion : Cell types per cluster                                 ")
message("+-------------------------------------------------------------------------------")

tmp_df3 <- as.data.frame(table(tmp_df[,c("CLUSTER", "Experiment")]))
tmp_df3 <- subset(tmp_df3, tmp_df3$Freq > 30)
tmp_df3 <- as.data.frame(table(tmp_df3$CLUSTER))
mean(tmp_df3$Freq)

cells_in_clusters <- as.data.frame(table(tmp_df$CLUSTER))


tmp_df4 <- as.data.frame(table(tmp_df[,c("CLUSTER", "Experiment")]))
tmp_df4$total_cells <- cells_in_clusters[match(tmp_df4$CLUSTER, cells_in_clusters$Var1),]$Freq
tmp_df4$perc_per_clust <- tmp_df4$Freq / tmp_df4$total_cells * 100
tmp_df4 <- subset(tmp_df4, tmp_df4$perc_per_clust > 5)
tmp_df4 <- as.data.frame(table(tmp_df4$CLUSTER))
mean(tmp_df4$Freq)

tmp_df4 <- as.data.frame(table(tmp_df[,c("CLUSTER", "Experiment")]))
tmp_df4$total_cells <- cells_in_clusters[match(tmp_df4$CLUSTER, cells_in_clusters$Var1),]$Freq
tmp_df4$perc_per_clust <- tmp_df4$Freq / tmp_df4$total_cells * 100
tmp_df4 <- subset(tmp_df4, tmp_df4$perc_per_clust > 10)
tmp_df4 <- as.data.frame(table(tmp_df4$CLUSTER))
mean(tmp_df4$Freq)

message("--------------------------------------------------------------------------------")
message("+           Dispersion : clusters per celltype                                  ")
message("+-------------------------------------------------------------------------------")

tmp_df5 <- as.data.frame(table(tmp_df[,c( "Experiment","CLUSTER")]))
tmp_df5 <- subset(tmp_df5, tmp_df5$Freq > 100)
tmp_df5 <- as.data.frame(table(tmp_df5$Experiment))
mean(tmp_df5$Freq)

length(unique(tmp_df$CLUSTER))

message("--------------------------------------------------------------------------------")
message("+          How many clusters have just 1 batch (>90% cells)                     ")
message("+-------------------------------------------------------------------------------")

tmp_df6 <- as.data.frame(table(tmp_df[,c( "CLUSTER", "Batch")]))
tmp_df6$total_cells <- cells_in_clusters[match(tmp_df6$CLUSTER, cells_in_clusters$Var1),]$Freq
tmp_df6$perc_per_clust <- tmp_df6$Freq / tmp_df6$total_cells * 100

tmp_df6_90 <- subset(tmp_df6, tmp_df6$perc_per_clust > 90)
tmp_df6_90 <- as.data.frame(table(tmp_df6_90$CLUSTER))
nrow(tmp_df6_90[tmp_df6_90$Freq == 1,])
nrow(tmp_df6_90)

tmp_df6_95 <- subset(tmp_df6, tmp_df6$perc_per_clust > 95)
tmp_df6_95 <- as.data.frame(table(tmp_df6_95$CLUSTER))
nrow(tmp_df6_95[tmp_df6_95$Freq == 1,])
nrow(tmp_df6_95)

tmp_df6_98 <- subset(tmp_df6, tmp_df6$perc_per_clust > 98)
tmp_df6_98 <- as.data.frame(table(tmp_df6_98$CLUSTER))
nrow(tmp_df6_98[tmp_df6_98$Freq == 1,])
nrow(tmp_df6_98)







