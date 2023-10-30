#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()
options(bitmapType='cairo')

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/TEMPORA"
setwd(baseDir)

library(Tempora)

dataset <- "R"
#dataset <- "I"
#dataset <- "Together"




message("+-------------------------------------------------------------------------------")
message("+                             Downsampling                                      ")
message("+-------------------------------------------------------------------------------")
library(dplyr)

#matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_clust24_25rm.Rds")
matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_matrix.su_R__29__umap_n.n25_repuls2_spr2L.Rds")
#matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_matrix.su_I__32__umap_n.n25_repuls2_spr2L.Rds")

gc()
nrow(matrix.su@meta.data)

table(matrix.su@meta.data[,c("Experiment", "Batch")])
matrix.su@meta.data$Experiment_batch <- paste(matrix.su@meta.data$Experiment, matrix.su@meta.data$Batch, sep = "_")

tmp_tbl <- as.data.frame(matrix.su@meta.data$Experiment_batch)
rownames(tmp_tbl) <- rownames(matrix.su@meta.data)
colnames(tmp_tbl)[1] <- "Experiment"
tbl_order <- as.data.frame(table(tmp_tbl))
tbl_order <- tbl_order[order(tbl_order$Freq),]
Experiment_batch_names <-  as.character(tbl_order$tmp_tbl)

cells_to_keep <- list()
for (i in 1:length(Experiment_batch_names)) {
  sample_expt_tmp <- subset(tmp_tbl, tmp_tbl$Experiment == Experiment_batch_names[i])
  sample_expt_tmp$cells <- rownames(sample_expt_tmp)
  if (nrow(sample_expt_tmp) > 5000){
    cells_to_keepp <- sample_n(sample_expt_tmp,  5000 ) #  as.integer(nrow(sample_expt_tmp)*0.8)
    cells_to_keep[[i]] <- cells_to_keepp$cells
  } else {  cells_to_keep[[i]] <- sample_expt_tmp$cells  }
}
cells_to_keep <- unlist(cells_to_keep)
length(cells_to_keep)

matrix.su.downsample <- subset(matrix.su, cells = cells_to_keep)

#saveRDS(matrix.su.downsample, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/DOWNSAMPLED_DATASETS/harmony_matrix.su_R__29__umap_n.n25_repuls2_spr2L_DOWNSAMPLED_35k.Rds")
#saveRDS(matrix.su.downsample, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/DOWNSAMPLED_DATASETS/harmony_matrix.su_I__32__umap_n.n25_repuls2_spr2L_DOWNSAMPLED_42k.Rds")
#saveRDS(matrix.su.downsample, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/DOWNSAMPLED_DATASETS/harmony_matrix.su_together__umap_n.n25_repuls2_spr2L_DOWNSAMPLED_72k.Rds")






message("+-------------------------------------------------------------------------------")
message("+                       read in dataset                                         ")
message("+-------------------------------------------------------------------------------")

dataset <- "R"

if(dataset == "I"){
  matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/DOWNSAMPLED_DATASETS/harmony_matrix.su_I__32__umap_n.n25_repuls2_spr2L_DOWNSAMPLED_42k.Rds")
  print("dataset:: I")
  Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.0.8
  matrix.su@meta.data$Cluster <- matrix.su@meta.data$SCT_snn_res.0.8
} else if (dataset == "R") {
  matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/DOWNSAMPLED_DATASETS/harmony_matrix.su_R__29__umap_n.n25_repuls2_spr2L_DOWNSAMPLED_35k.Rds")
  print("dataset:: R")
  Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.0.8
  matrix.su@meta.data$Cluster <- matrix.su@meta.data$SCT_snn_res.0.8
} else if (dataset == "Together"){
  matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/DOWNSAMPLED_DATASETS/harmony_matrix.su_together__umap_n.n25_repuls2_spr2L_DOWNSAMPLED_72k.Rds")
  print("dataset:: Together")
  Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.1
  matrix.su@meta.data$Cluster <- matrix.su@meta.data$SCT_snn_res.1
}else{ print("no dataset assigned") }

head(matrix.su@meta.data,3)

matrix.su@meta.data$Age <- matrix.su@meta.data$Experiment
matrix.su@meta.data$Age <- gsub(  "*.\\.", "", matrix.su@meta.data$Age )
matrix.su@meta.data$Age <- factor(matrix.su@meta.data$Age, levels = c("0", "1", "4", "24", "36", "48"))

unique(matrix.su@meta.data$Cluster)
table(matrix.su@meta.data$Cluster)

head(matrix.su@meta.data,2)
colnames(matrix.su@meta.data)

head(matrix.su@meta.data,2)

is.numeric(matrix.su@meta.data$cluster_tempora)

matrix.su@meta.data$cluster_tempora <- as.numeric(as.character(matrix.su@meta.data$Cluster))
Idents(matrix.su) <- matrix.su@meta.data$cluster_tempora

table(matrix.su@meta.data$cluster_tempora)

cluster_string <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23", "24", "25")

matrix.su@meta.data$cluster_tempora <- factor(matrix.su@meta.data$cluster_tempora, levels = cluster_string[1:length(unique(matrix.su@meta.data$cluster_tempora))])


#  Import matrix.su_I_downsampled to 42000c ::: 
tempora_obj <- ImportSeuratObject(matrix.su, assayType = "SCT",  clusters = "cluster_tempora",
                                  timepoints = "Age", 
                                  cluster_labels = levels(matrix.su@meta.data$cluster_tempora),
                                  timepoint_order = c("0",  "1", "4" ,"24", "36",  "48")) # 

head(tempora_obj@meta.data,2)
head(tempora_obj@cluster.metadata,3)
tempora_obj@cluster.metadata <- tempora_obj@cluster.metadata[order(tempora_obj@cluster.metadata$Cluster_time_score),]

tempora_obj@data[1:5,1:5]
dim(tempora_obj@data) 
dim(tempora_obj@meta.data) 
tempora_obj@cluster.metadata$Id

#Estimate pathway enrichment profiles of clusters
#tempora_obj <- CalculatePWProfiles(tempora_obj, gmt_path = "Mouse_GOBP_AllPathways_no_GO_iea_February_01_2020_symbol.gmt", method="ssgsea", min.sz = 5, max.sz = 200, parallel.sz = 1)
tempora_obj <- CalculatePWProfiles(tempora_obj, gmt_path = "Mouse_WikiPathways_February_01_2020_symbol.gmt", method="ssgsea", min.sz = 5, max.sz = 500, parallel.sz = 1) # dont use gsva  Gene Set Variation Analysis (GSVA)1    ssgsea

head(tempora_obj@cluster.pathways)
sum(rownames(tempora_obj@meta.data) == colnames(tempora_obj@data)) # 24763

#Build trajectory with 6 PCs 
#We can now build the trajectory based on the clusters’ pathway enrichment profiles. Tempora employs the mutual information (MI) rank and data processing inequality approach implemented in ARACNE to calculate MI between all cluster pairs present in the data as well as remove edges with weak MIs. The trajectory is stored as a dataframe of edge lists in the trajectory slot. Tempora then assigns directions to all edges in the network so that edges point from clusters with low temporal scores to clusters with high temporal scores.
tempora_obj <- BuildTrajectory(tempora_obj, n_pcs = 6, difference_threshold = 0.01)


#Visualize the trajectory
tempora_obj <- PlotTrajectory(tempora_obj)
#Error in igraph::graph_from_data_frame(d = object@trajectory, vertices = object@cluster.metadata,  : Some vertex names in edge list are not listed in vertex data frame
#tempora_obj@cluster.metadata$check <- tempora_obj@cluster.metadata$label
#tempora_obj@cluster.metadata$check <- gsub("Cluster " ,"" , tempora_obj@cluster.metadata$check)
#head(tempora_obj@cluster.metadata)
#tempora_obj@trajectory$check <- paste(tempora_obj@trajectory$from, tempora_obj@trajectory$to, sep = "-")
tempora_obj@trajectory
tempora_obj@cluster.metadata


library(dplyr)
tempora_obj@trajectory[which(! tempora_obj@trajectory$from %in% tempora_obj@cluster.metadata$Id) ,]
tempora_obj@trajectory[which(! tempora_obj@trajectory$to %in% tempora_obj@cluster.metadata$Id) ,]

tempora_obj@trajectory <- tempora_obj@trajectory[-which(! tempora_obj@trajectory$from %in% tempora_obj@cluster.metadata$Id) ,]
tempora_obj@trajectory <- tempora_obj@trajectory[-which(! tempora_obj@trajectory$to %in% tempora_obj@cluster.metadata$Id) ,]

dim(tempora_obj@trajectory)

tempora_obj <- PlotTrajectory(tempora_obj, curved = TRUE, vertex.label.dist = 1)


pdf(paste( "TEMPORA_trajectory", dataset, "downsampled", ".pdf", sep="_"), width=16,height=12)
par(bg=NA)
PlotTrajectory(tempora_obj, curved = TRUE, vertex.label.dist = 1, label.degree = 45, label.font= 'sans')
dev.off()



#tempora_obj@cluster.metadata <- tempora_obj@cluster.metadata[tempora_obj@cluster.metadata$check %in% tempora_obj@trajectory$check,]
#dim(tempora_obj@cluster.metadata)
#dim(tempora_obj2@cluster.metadata)


#FIND TIME DEPENDENT PATHWAYS ::: Fit GAMs on pathway enrichment profile
tempora_obj <- IdentifyVaryingPWs(tempora_obj, pval_threshold = 0.05)

tempora_obj@trajectory
tempora_obj@cluster.metadata$Id


#  R :::: 16 18 13 1  2  20 6  14 0  3  15 21 11 5  9  4  17 10 8  7  12 19
#  I :::: 9  6  1  17 14 8  11 0  7  3  4  12 13 10 18 15 2  5  16
# together:::::   12 6  18 1  19 10 23 16 2  14 11 5  15 0  9  17 8  13 4  20 3  21 7  22