#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()
options(bitmapType='cairo')


library(cowplot)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 4000 * 1024^2)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
Project <- "CTR_dscj1_HARMONY_markers"
#Project <- "CTR_dscj1_IntegrationSCT_wRef"




matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")

nrow(matrix.su@meta.data)             # 229793 
min(matrix.su@meta.data$nFeature_RNA) # 300
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N721",])
matrix.su <- subset(matrix.su, orig.ident2 != "Batch0004_9617_N721")
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N729")

head(matrix.su@meta.data,2)
tmp_df <- unique(matrix.su@meta.data[,c("Experiment", "orig.ident2", "Batch")])
table(tmp_df[,c(1,3)])






message("+-------------------------------------------------------------------------------")
message("+                                  HARMONY                                       ")
message("+-------------------------------------------------------------------------------")

#matrix.su <- NormalizeData(matrix.su)
#matrix.su <- FindVariableFeatures(matrix.su)
#matrix.su <- ScaleData(matrix.su)
#matrix.su <- RunPCA(matrix.su, npcs = 20, verbose = FALSE)

harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT")
#harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", theta = 4, reference_values = "0B0")

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01)
harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony", dims = 1:30) 
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1.5)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.6)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.8)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0)

#saveRDS(harmony_matrix.su, "harmony_theta4_ref0B0_matrix.su.Rds")

#saveRDS(harmony_matrix.su, "harmony_matrix.su.Rds")

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 20)
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_n.neigh20.Rds")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 40)
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_n.neigh40.Rds")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, local.connectivity = 5L)
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_local.conn5L.Rds")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, spread = 3)
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_spread3L.Rds")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 25)
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_n.neigh25.Rds")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 50)
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_n.neigh50.Rds")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 10)
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_n.neigh10.Rds")
harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 25, repulsion.strength = 2, spread = 2L )
#saveRDS(harmony_matrix.su, "harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729.Rds")


##dropseq.integrated <- readRDS("harmony_matrix.su.Rds")
##dropseq.integrated <- readRDS("harmony_matrix.su_umap_n.neigh10.Rds")
#dropseq.integrated <- readRDS("harmony_matrix.su_umap_n.neigh20.Rds")
#dropseq.integrated <- readRDS("harmony_matrix.su_umap_n.neigh25.Rds")
##dropseq.integrated <- readRDS("harmony_matrix.su_umap_n.neigh40.Rds")
#dropseq.integrated <- readRDS("harmony_matrix.su_umap_n.neigh50.Rds")
#dropseq.integrated <- readRDS("harmony_matrix.su_umap_local.conn5L.Rds")
#dropseq.integrated <- readRDS("harmony_matrix.su_umap_spread3L.Rds")
#dropseq.integrated <- readRDS("HARMONY/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L.Rds")
dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729.Rds")
gc()

#dropseq.integrated <- readRDS("HARMONY/harmony_theta4_ref0B0_matrix.su.Rds")
#dropseq.integrated <- readRDS("HARMONY/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L.Rds")
#dropseq.integrated <- FindNeighbors(dropseq.integrated, reduction = "harmony", dims = 1:30) 
#dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 1.25)



message("+-------------------------------------------------------------------------------")
message("+                      Seurat Integration  with REFERENCE                       ")
message("+-------------------------------------------------------------------------------")
dropseq.features <- readRDS("dropseq.features.Rds")
dropseq.list <- readRDS("dropseq.list_prepped.Rds")

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(dropseq.list) == "0B0")

dropseq.anchors <- FindIntegrationAnchors(object.list = dropseq.list, normalization.method = "SCT", anchor.features = dropseq.features, reference = reference_dataset)

dropseq.integrated <- IntegrateData(anchorset = dropseq.anchors, normalization.method = "SCT")
gc()
dropseq.integrated <- RunPCA(dropseq.integrated, verbose = TRUE)
gc()
dropseq.integrated <- RunUMAP(dropseq.integrated, dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L)
gc()
#saveRDS(dropseq.integrated, "dropseq.integrated_SCT_Ref0B0_clust.Rds")
gc()

dropseq.integrated <- readRDS("dropseq.integrated_SCT_Ref0B0_clust.Rds")




message("+-------------------------------------------------------------------------------")
message("+                      Seurat Integration                                       ")
message("+-------------------------------------------------------------------------------")

dropseq.list <- SplitObject(matrix.su, split.by = "Batch")
names(dropseq.list) # "ABC" "0B0" "AB0" "A0C" "0BC" "A00" "00C" ---> 0B0 to use as reference as has t0 and more timepoints than others.
for (i in 1:length(dropseq.list)) {
  dropseq.list[[i]] <- SCTransform(dropseq.list[[i]], verbose = TRUE)
}

#saveRDS(dropseq.list, "dropseq.list.Rds")

#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.

dropseq.features <- SelectIntegrationFeatures(object.list = dropseq.list, nfeatures = 3000)
#saveRDS(dropseq.features, "dropseq.features.Rds")

dropseq.list <- PrepSCTIntegration(object.list = dropseq.list, anchor.features = dropseq.features, verbose = TRUE)
#gc()
#saveRDS(dropseq.list, "dropseq.list_prepped.Rds")



#Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set normalization.method = 'SCT':
  
dropseq.anchors <- FindIntegrationAnchors(object.list = dropseq.list, normalization.method = "SCT", anchor.features = dropseq.features, verbose = TRUE)
#saveRDS(dropseq.anchors, "dropseq.anchors.Rds")

dropseq.integrated <- IntegrateData(anchorset = dropseq.anchors, normalization.method = "SCT",  verbose = TRUE)
#saveRDS(dropseq.integrated, "dropseq.integrated.Rds")


#Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset. Commands are identical to the standard workflow, but do not run the ScaleData function after integration. You can see that after integration, cells group by their biological cell type (which has been pre-annotated), instead of by their underlying technology.

dropseq.integrated <- RunPCA(dropseq.integrated, verbose = TRUE)
#dropseq.integrated <- RunUMAP(dropseq.integrated, dims = 1:30)
dropseq.integrated <- RunUMAP(dropseq.integrated, dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 


dropseq.integrated <- FindNeighbors(dropseq.integrated, dims = 1:30)
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 0,    dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1)  
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 0.2,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 0.6,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 0.8,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 1,    dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 1.25,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 1.5,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
dropseq.integrated <- FindClusters(dropseq.integrated, resolution = 2,    dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 


#
#saveRDS(dropseq.integrated, "dropseq.integrated_clust.Rds")
dropseq.integrated <- readRDS("SCT_integration/dropseq.integrated_clust_CC.Rds")
colnames(dropseq.integrated@meta.data)

message("+-------------------------------------------------------------------------------")
message("+                   Add missing metedata!!!!                                    ")
message("+-------------------------------------------------------------------------------")

colnames(dropseq.integrated@meta.data)
#colnames(matrix.su@meta.data)
#phase_assignment <- matrix.su@meta.data[,c(1,2,5,9,21,22,23)]
#saveRDS(phase_assignment, "CellCycle_phase_assignment.Rds")

phase_assignment <- readRDS("CellCycle_phase_assignment.Rds")
dropseq.integrated <- readRDS("dropseq.integrated_clust.Rds")

dropseq.integrated@meta.data$Phase      <- phase_assignment[match(rownames(dropseq.integrated@meta.data), rownames(phase_assignment)),]$Phase
dropseq.integrated@meta.data$G2M.Score  <- phase_assignment[match(rownames(dropseq.integrated@meta.data), rownames(phase_assignment)),]$G2M.Score
dropseq.integrated@meta.data$S.Score    <- phase_assignment[match(rownames(dropseq.integrated@meta.data), rownames(phase_assignment)),]$S.Score
dropseq.integrated@meta.data$percent.mt <- phase_assignment[match(rownames(dropseq.integrated@meta.data), rownames(phase_assignment)),]$percent.mt

#saveRDS(dropseq.integrated, "dropseq.integrated_clust_CC.Rds")

  

#dropseq.integrated <- readRDS("dropseq.integrated_clust_CC.Rds")

  



message("+-------------------------------------------------------------------------------")
message("+              Explore and QC integrated dataset                                ")
message("+-------------------------------------------------------------------------------")
  
plots <- DimPlot(dropseq.integrated, group.by = c("Batch", "Experiment"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))

plots2 <- DimPlot(dropseq.integrated, group.by = "Batch", split.by = "Experiment", ncol = 3)



message("+-------------------------------------------------------------------------------")
message("+                       CLUSTREE                                                ")
message("+-------------------------------------------------------------------------------")


library(clustree)
tree_to_plot <- clustree::clustree(dropseq.integrated, prefix = "integrated_snn_res.") # "RNA_snn_res."
tree_to_plot <- clustree::clustree(dropseq.integrated, prefix = "SCT_snn_res.") # "RNA_snn_res."

#pdf(paste("dropseq.integrated_clustree.pdf", sep=""), width=12,height=8)
pdf(paste("harmony_matrix.su_clustree.pdf", sep=""), width=12,height=8)
par(bg=NA)
tree_to_plot
dev.off()






message("+-------------------------------------------------------------------------------")
message("+                        Plotting UMAP - custom                                 ")
message("+-------------------------------------------------------------------------------")

sampleTable <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/sampleTable_MERGED_53.csv")
dataset_name <- "Harmony_25_rep2_spr2_n729rm" #"dropseq.integrated"
head(dropseq.integrated@meta.data,2)

matrix.umap            <- as.data.frame(Embeddings(object=dropseq.integrated, reduction="umap"))
matrix.umap$Sample_2   <- rownames(matrix.umap) 
matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
rownames(sampleTable)  <- sampleTable$Sample_2 
matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
matrix.umap$Project    <- sampleTable[matrix.umap$Sample_2, ]$Project
head(matrix.umap)

cell_cycle             <- dropseq.integrated@meta.data %>% dplyr::select((Phase))  # 
colnames(cell_cycle)   <- c("Phase")
matrix.umap$Phase      <- cell_cycle$Phase

#reso                <- "SCT_snn_res.0.8"
#clust               <- dropseq.integrated@meta.data %>% dplyr::select((integrated_snn_res.1))  # 
reso                <- "SCT_snn_res.1.0"
clust               <- dropseq.integrated@meta.data %>% dplyr::select((SCT_snn_res.1))  # 

colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 

head(matrix.umap,2)
unique(matrix.umap$Batch)
unique(matrix.umap$Cluster)


table(matrix.umap$Experiment)
table(matrix.umap$Cluster)

message("----------------------Cluster annotation -----------------------")

Clusters <- as.numeric(levels(matrix.umap$Cluster))

coord_x <- list()
coord_y <- list()

for (i in 1:length(Clusters)){
  coord_x[i] <- median(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$UMAP_1)  
  coord_y[i] <- median(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$UMAP_2)  
}

head(coord_x)
names(coord_x) <- Clusters
names(coord_y) <- Clusters

Cluster_Anno <- list()

for (i in 1:length(Clusters) ){
  #Cluster_Anno[i] <-  paste0(  "annotate(\'text\', label = paste0(\'Cluster_\'," , Clusters[i], "), x =", coord_x[i], ", y = ", coord_y[i], ", color = \'black\')" )
  #Cluster_Anno[i] <-  paste0(  "annotate(\'text\', label = paste0(", Clusters[i], "), x =", coord_x[i], ", y = ", coord_y[i], ", color = \'black\', size=8)" )
  Cluster_Anno[i] <-  paste0(  "annotate(\'text\', label = paste0(\'", Clusters[i], "\'), x =", coord_x[i], ", y = ", coord_y[i], ", color = \'black\', size=8)" )
  
}

Cluster_Anno <- unlist(Cluster_Anno)

library(stringr)
Cluster_Anno <- str_c(Cluster_Anno,  collapse =  " + ")






message("----------------------umap.Cluster-----------------------")
unique(matrix.umap$Cluster)
matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22","23", "24", "25")) # "23", "24", "25", "26" ))
#matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29" ))

umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.2, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0(dataset_name, " QC :::  ", reso, " Cluster" )) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() +
  annotate('text', label = paste0('0'), x =-11.1303195138008, y = 0.448520670055307, color = 'black', size=8) + annotate('text', label = paste0('1'), x =5.21850570249256, y = 0.588521832941927, color = 'black', size=8) + annotate('text', label = paste0('10'), x =1.33514096545872, y = 9.54303281974689, color = 'black', size=8) + annotate('text', label = paste0('11'), x =14.1224566321343, y = 3.8776468582143, color = 'black', size=8) + annotate('text', label = paste0('12'), x =12.0221615653008, y = -2.82656246948346, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-12.1390885491401, y = -5.06648117828473, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.62632545041737, y = 2.00741069984332, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-13.9796017785102, y = 4.23234718513385, color = 'black', size=8) + annotate('text', label = paste0('16'), x =10.188165746209, y = 8.16320484352008, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0143977026909, y = 14.322745973586, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-11.977231897834, y = 1.29518311691181, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77784530115429, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.29146710825268, y = -4.13745767402752, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.69010885762867, y = 3.4951583690633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('22'), x =2.96324016856846, y = -0.739922022105299, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.0985351424187, y = 4.13216405582324, color = 'black', size=8) + annotate('text', label = paste0('4'), x =7.88796337651905, y = -0.811322336674772, color = 'black', size=8) + annotate('text', label = paste0('5'), x =8.10612925100025, y = -8.10546809959515, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-0.444073416712914, y = -7.71475559997662, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-11.1568364281685, y = -0.251747911454283, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.08930643605884, y = 6.20413273048297, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.884592451575433, y = 0.91303431463138, color = 'black', size=8)

umap.Cluster

#umap.Cluster_1 <- umap.Cluster
#umap.Cluster_0.8 <- umap.Cluster
#umap.Cluster_0.6 <- umap.Cluster

#umap.Cluster_0.3 <- umap.Cluster
#umap.Cluster_0.4 <- umap.Cluster
#umap.Cluster_0.5 <- umap.Cluster

png(paste0("UMAP30_MergedRealigned_51_ini_", Project,  "_ggplot_", "Cluster", reso,".png"), width=1000, height=900, type="cairo")
#pdf(paste0("UMAP30_MergedRealigned_53_ini_", Project, "_ggplot_", "Expt", ".pdf"), width=10,height=9)
par(bg=NA)
umap.Cluster
dev.off()





message("----------------------umap.QC_1_sample-----------------------")

table(dropseq.integrated@meta.data[,c(5,8)])

matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_samples")
matrix.umap$BatchSample <- paste( matrix.umap$Batch, matrix.umap$Sample, sep = "_")
unique(matrix.umap[matrix.umap$QC == "0.0",]$BatchSample)
# "0B0_Batch0011_15988_N704" "0BC_Dups_setA1_reps_N724" "A00_Dups_setA2_reps_N719" "A00_Dups_setB_reps_N704"  "00C_Batch0004_9617_N722" 
# "00C_Batch0004_9617_N726"  "00C_Batch0004_9617_N727"  "00C_Batch0004_9617_N720"  "00C_Batch0004_9617_N719"  "00C_Batch0004_9617_N718" 
#  "00C_Batch0004_9617_N729" 

matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "0B0_Batch0011_15988_N704", "0B0_Batch0011_15988_N704", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "0BC_Dups_setA1_reps_N724", "0BC_Dups_setA1_reps_N724", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "A00_Dups_setA2_reps_N719", "A00_Dups_setA2_reps_N719", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "A00_Dups_setB_reps_N704", "A00_Dups_setB_reps_N704", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N722", "00C_Batch0004_9617_N722", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N726", "00C_Batch0004_9617_N726", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N727", "00C_Batch0004_9617_N727", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N720", "00C_Batch0004_9617_N720", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N719", "00C_Batch0004_9617_N719", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N718", "00C_Batch0004_9617_N718", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N729", "00C_Batch0004_9617_N729", matrix.umap$QC)




#matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "A00", "A00_t0", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "00C", "00C_t0", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0B0", "0B0_t0", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0BC", "0BC_t0", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.3, size=0.1) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  #scale_colour_manual("", values = c("other_samples"="lightgrey", "A00_t0" = "blue1", "00C_t0"="red", "0B0_t0" = "orange", "0BC_t0" = "green")) +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "0B0_Batch0011_15988_N704" = "blue1", "0BC_Dups_setA1_reps_N724"="red", "A00_Dups_setA2_reps_N719" = "orange", "A00_Dups_setB_reps_N704" = "green", "00C_Batch0004_9617_N722"= "khaki3", "00C_Batch0004_9617_N726" = "brown", "00C_Batch0004_9617_N727" = "navyblue", "00C_Batch0004_9617_N720"= "pink", "00C_Batch0004_9617_N719"= "violet", "00C_Batch0004_9617_N718"= "turquoise")) + # , "00C_Batch0004_9617_N729" = "black"
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() #+ xlim(c(6,7)) + ylim(c(-3, -1.8)) # + xlim(c(3,9)) + ylim(c(-3.5, 2.5))
umap.qc


matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1", "I.1", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "0B0", "0B0_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "A00", "A00_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "A0C", "A0C_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "ABC", "ABC_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "AB0", "AB0_I.1", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.1 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "A00_I.1" = "blue1", "A0C_I.1"="red", "0B0_I.1" = "orange", "ABC_I.1" = "green", "AB0_I.1"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24", "I.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$Batch == "ABC", "ABC_I.24", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$Batch == "AB0", "AB0_I.24", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.24 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_I.24" = "green", "AB0_I.24"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc



matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4", "R.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "ABC", "ABC_R.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "AB0", "AB0_R.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "A0C", "A0C_R.4", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC r.4 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_R.4" = "green", "AB0_R.4"="purple", "A0C_R.4"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4", "I.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "ABC", "ABC_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "AB0", "AB0_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "A0C", "A0C_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "A00", "A00_I.4", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC i.4 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_I.4" = "green", "AB0_I.4"="purple", "A0C_I.4"="blue", "A00_I.4" = "orange")) +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24", "R.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24" & matrix.umap$Batch == "ABC", "ABC_R.24", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24" & matrix.umap$Batch == "A00", "A00_R.24", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC r.24 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_R.24" = "green", "A00_R.24"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc




#png(paste0("UMAP30_MergedRealigned_53_ini_", Project, "_ggplot_", "CellCycle", ".png"), width=1000,height=900, type="cairo")
#par(bg=NA)
#umap.cc
#dev.off()




message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " Cell Cycle" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() #   theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +

umap.cc

png(paste0("UMAP_dropseq", Project, "theta4_ref0B0", "_ggplot_", "CellCycle", ".png"), width=1000,height=900, type="cairo")
par(bg=NA)
umap.cc
dev.off()


message("----------------------umap.Experiment-----------------------")

umap.Experiment    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.2, size=0.3) +
  ggtitle(paste0("Resolution ", reso, " :: Experiment " )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Experiment, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="black", "I.1"="yellow", "R.1"="yellow2", 
                                     "I.4"="orangered", "R.4"="red", "I.24"="lightblue2", "R.24"="lightblue4",
                                     "I.36"="blue", "R.36"="darkblue","I.48"="purple", "R.48"="orchid"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  theme_classic() + coord_fixed() +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=6), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-11.2057646889717, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.84784122514423, y = 2.00883030366794, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-12.036913790229, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + annotate('text', label = paste0('24'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('25'), x =3.5460153441399, y = -0.737795358182036, color = 'black', size=8)

#
umap.Experiment

png(paste0("UMAP_MergedRealigned_51_int_", Project, "_ggplot_", "Expt", ".png"), width=1000, height=900, type="cairo")
#pdf(paste0("UMAP30_MergedRealigned_53_ini_", Project, "_ggplot_", "Expt", ".pdf"), width=10,height=9)
par(bg=NA)
umap.Experiment
dev.off()








message("----------------------umap.t.0-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_timepoints")
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "A00", "A00_t0", matrix.umap$Mixer)
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "00C", "00C_t0", matrix.umap$Mixer)
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0B0", "0B0_t0", matrix.umap$Mixer)
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0BC", "0BC_t0", matrix.umap$Mixer)
table(matrix.umap$Mixer)

umap.t0    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.2, size=0.2) +
  ggtitle(paste0("Resolution ", reso, " T0 in Batches" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "A00_t0" = "blue1", "00C_t0"="red", "0B0_t0" = "orange", "0BC_t0" = "green")) +
  coord_fixed() +
  xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-11.2057646889717, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.84784122514423, y = 2.00883030366794, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-12.036913790229, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + annotate('text', label = paste0('24'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('25'), x =3.5460153441399, y = -0.737795358182036, color = 'black', size=8)




# clusters 15,23,24
head(dropseq.integrated@meta.data,2)
#table(dropseq.integrated@meta.data[dropseq.integrated@meta.data$Experiment == "0.0" , c("Experiment", "SCT_snn_res.1")])
table(dropseq.integrated@meta.data[ , c("Experiment", "SCT_snn_res.1")])
table(dropseq.integrated@meta.data[dropseq.integrated@meta.data$Experiment == "0.0" , c("Batch", "SCT_snn_res.1")])



message("----------------------umap.t.1-----------------------")

umap.Experiment1    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="blue", "R.1"="red", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-11.2057646889717, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.84784122514423, y = 2.00883030366794, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-12.036913790229, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + annotate('text', label = paste0('24'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('25'), x =3.5460153441399, y = -0.737795358182036, color = 'black', size=8)




message("----------------------umap.t.4-----------------------")

umap.Experiment4    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="blue", "R.4"="red", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-11.2057646889717, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.84784122514423, y = 2.00883030366794, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-12.036913790229, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + annotate('text', label = paste0('24'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('25'), x =3.5460153441399, y = -0.737795358182036, color = 'black', size=8)




message("----------------------umap.t.24-----------------------")

umap.Experiment24    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="blue", "R.24"="red",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c('0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))  + 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-11.2057646889717, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.84784122514423, y = 2.00883030366794, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-12.036913790229, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + annotate('text', label = paste0('24'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('25'), x =3.5460153441399, y = -0.737795358182036, color = 'black', size=8)






message("----------------------umap.t.36-----------------------")

umap.Experiment36    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="blue", "R.36"="red","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-11.2057646889717, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.84784122514423, y = 2.00883030366794, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-12.036913790229, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + annotate('text', label = paste0('24'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('25'), x =3.5460153441399, y = -0.737795358182036, color = 'black', size=8)







message("----------------------umap.t.48-----------------------")

umap.Experiment48    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="blue", "R.48"="red"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-11.2057646889717, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.84784122514423, y = 2.00883030366794, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-12.036913790229, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + annotate('text', label = paste0('24'), x =1.15898526417908, y = -4.07302409935101, color = 'black', size=8) + annotate('text', label = paste0('25'), x =3.5460153441399, y = -0.737795358182036, color = 'black', size=8)



png(paste( "UMAP_MergedRealigned_51_int_", Project, "rmN729","ggplot", "Expt_grid_eachTimePoint", reso, ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
plot_grid(umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 2, align = "hv")
dev.off()








message("+-------------------------------------------------------------------------------")
message("+                 plot previously identified genes                              ")
message("+-------------------------------------------------------------------------------")

# top QC check genes: 
FeaturePlot(dropseq.integrated, features = c("Elf5", "Eomes", "Sox2","Ascl2","Gcm1","Ovol2", "Anxa1" , "Plac1", "Uba6")) # "Anxa1,"Pparg" "Cdx2"

png(paste( "FeaturePlt_UMAP_MergedRealigned_52_int_", Project, "theta4_ref0B0","ggplot",  ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
FeaturePlot(dropseq.integrated, features = c("Elf5", "Eomes", "Sox2","Ascl2","Gcm1","Ovol2", "Anxa1" , "Plac1", "Uba6")) # "Anxa1,"Pparg" "Cdx2"
dev.off()



FeaturePlot(dropseq.integrated, features = c("Elf5", "Cdx2","Eomes","Tead4","Sox2","Ascl2","Gcm1","Ovol2"))
# LP ("Gcm1","Ovol2"), EPC = Ascl2


FeaturePlot(dropseq.integrated, features = c("Usf2", "Nfya", "Yy1", "Maff", "Gcm1","Ovol2", "Pparg","Ghrl1")) # LP
FeaturePlot(dropseq.integrated, features = c("Vezf1", "Gata1", "Klf7", "Bdp1","Ascl2")) # EPC
FeaturePlot(dropseq.integrated, features = c("Bhlhe40", "Smarca4", "Myc", "Foxo3","Sp3","Srebf2","Hdac6","Tead3", "Ppard")) # TFs
# other TFs : "E2f8","","Grhl1","","Foxo4","Hdac2","","Junb","Nfkb2","Creb3l2","Pou2f1","Mycn","","Max","","","E2f7","Fos","Elf2"


VlnPlot(dropseq.integrated, features = c("Lgals3", "Krt18"), assay = "SCT", log = TRUE)
FeaturePlot(dropseq.integrated, features = c("Flt1", "Anxa1", "Uba6", "Bnip3","Tceb2", "Plac1")) # earlier markrs




# T0 identified genes: 
FeaturePlot(dropseq.integrated, features = c("Krt18", "Lgals3", "Malat1","Rhox6","Atrx","Cenpf", "Ubb" , "Nars", "Ncl")) # "Anxa1,"Pparg" "Cdx2"
FeaturePlot(dropseq.integrated, features = c("Rhox9", "Rpl37", "Krt8","Crip1","Igf2","Cdkn1c", "Npm1" , "Hand1", "Rpl41")) # "Anxa1,"Pparg" "Cdx2"






message("+-------------------------------------------------------------------------------")
message("+                         MARKERS for res 1.0                                   ")
message("+-------------------------------------------------------------------------------")

markerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/"
setwd(markerDir)

reso <- "res.1.0"
l2fc_cutoff <- 0.25
#result_list <- readRDS("HARMONY_together_MARKERS/result_list_harmony_together_res.1.0_l2fc0.6.Rds")
#marker_tbl <- readRDS( "HARMONY_together_MARKERS/marker_tbl.Rds")
result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/result_list_harmony_together_res.1.0_l2fc0.25.Rds")
marker_tbl <- readRDS( "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/marker_tbl_res.1.0_l2fc0.25.Rds")
last_cluster <- 25


get_res_tables <- function(x,y){
  if (is.data.frame(y[[x]]$result) == FALSE){ print(0)
  } else {
    y[[x]] <- y[[x]]$result
    tmp_df <- y[[x]]
    tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
    tmp_df <- subset(tmp_df, abs(tmp_df$avg_logFC) > l2fc_cutoff)
    return(nrow(tmp_df)) 
  }
}

calc_vals_for_marker_tbl <- function(x){
  marker_tbl_vals <- list()
  for(i in 1:length(x)){
    marker_tbl_vals[i] <- get_res_tables(i,x)
  }
  return(unlist(marker_tbl_vals))
}


marker_tbl_list <- list()
for (z in seq_along(result_list)){
  marker_tbl_list[[z]] <- calc_vals_for_marker_tbl(result_list[[z]])
  length(marker_tbl_list[[z]])
}
marker_tbl_list[[1]]


marker_tbl <- as.data.frame(marker_tbl)
length(marker_tbl_list)
for (z in seq_along(marker_tbl_list)){
  marker_tbl[z,] <- c( rep(NA, (last_cluster+1)-length(marker_tbl_list[[z]]) ), marker_tbl_list[[z]])
}
marker_tbl

length(marker_tbl_list[[1]])
length(marker_tbl_list[[2]])
length(marker_tbl_list[[3]])
length(marker_tbl_list[[4]])

max(marker_tbl, na.rm = T) # 78



library(biomaRt)
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene_id'), mart = ensembl)          



save_res_tables <- function(x,y){
  name_y <- res_name
  name_y <- gsub( "result_","" , name_y)
  name_x <- (x - 1 + as.numeric(name_y))
  if (is.data.frame(y[[x]]$result) == FALSE){ print(0)
  } else {
    y[[x]] <- y[[x]]$result
    if(y[[x]] == "same cluster"){  print(0)
    } else { 
      tmp_df <- y[[x]]
      tmp_df <- unique(merge(tmp_df, ensEMBL2id, by.x = "row.names", by.y = "external_gene_name", all.x = TRUE))
      tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
      tmp_df <- subset(tmp_df, abs(tmp_df$avg_logFC) > l2fc_cutoff)
      write.csv(tmp_df, paste0("CTR_dscj1_NewAlign51_HARMONY" , reso, "_l2fc", l2fc_cutoff,"_Markers_tbl__clusters_", name_y, ".vs.", name_x, ".csv"))
      print(nrow(y[[x]])) }
  }
}

save_res_tables(x=1, y=result_01)
save_res_tables(x=4, y=result_01)


for (z in seq_along(result_list)){
  res <- result_list[[z]]
  res_name <- names(result_list)[[z]]
  for (w in 1:length(res)){
    save_res_tables(w, res)
  }
}




library(ComplexHeatmap)
library(circlize)


f3 = colorRamp2( c(0, 1, 2, 3,5), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 

ht3 = Heatmap(as.matrix(marker_tbl),  col = f3, row_title = "", column_title = paste0("Markers between cluster pairs_", reso,  "_absl2fc", l2fc_cutoff), show_row_names = TRUE, heatmap_legend_param = list(title = "Number of markers", legend_height = unit(8, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left",  show_heatmap_legend = TRUE) #  width = unit(10, "cm"),
ht3


pdf(paste("Fig__Pairwise_MARKERS", Project, "ComplexHeatmap",  reso, "_l2fc",l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

#write.csv(marker_tbl, "CTR_dscj1_HARMONY_Together_51_MARKERS___Pairwise_cluster_marker_tbl__res.1_l2f0.25.csv", quote = FALSE)
#saveRDS(marker_tbl, "marker_tbl_res.1.0_l2fc0.25.Rds")






message("--------------------------------------------------------------------------------")
message("+                     GO  using enrichR                               ")
message("+-------------------------------------------------------------------------------")

library("enrichR")
enrichR_DB <- as.data.frame(listEnrichrDbs())

#l2fc_cutoff <- 0.6
l2fc_cutoff <- 0.25

#markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/HARMONY_together_MARKERS_GO_res1.0_l2fc0.6"
markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/HARMONY_together_MARKERS_GO_res1.0_l2fc0.25"
setwd(markerGODir)

save_GO_res <- function(x,y,db){
  db.term <- paste0(db, ".Term")
  name_y <- res_name
  name_y <- gsub( "result_","" , name_y)
  name_x <- (x - 1 + as.numeric(name_y))
  comparison_name <- paste0(name_y, ".vs.", name_x)
  print(comparison_name)
  if (is.data.frame(y[[x]]$result) == TRUE){ 
    y[[x]] <- y[[x]]$result
    tmp_df <- y[[x]]
    tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
    tmp_df <- subset(tmp_df, abs(tmp_df[,2]) > l2fc_cutoff)
    if( nrow(tmp_df) > 1){  
      enrichR_RES <- as.data.frame(enrichr(rownames(tmp_df), databases = db))
      enrichR_RES <- subset(enrichR_RES, enrichR_RES[,4] < 0.05)
      print(nrow(enrichR_RES))
      if( nrow(enrichR_RES) > 0 ){
        write.csv(enrichR_RES, paste0("NewAlign51_SCT_HARMONY__Res.1.0__", db, "_l2fc", l2fc_cutoff,"_", comparison_name, ".csv"))
        GO_matrix_to_add <- enrichR_RES[,c(1,4)]
        colnames(GO_matrix_to_add)[2] <- comparison_name
        GO_matrix <- unique(merge(GO_matrix, GO_matrix_to_add, by.x = db.term, by.y =db.term, all.x = TRUE, all.y = TRUE ))
        #print(head(GO_matrix,3))
        print(dim(GO_matrix))
      } 
    }
  }
  return(GO_matrix)
}


FeaturePlot(matrix.su, features = c("Saa3", "Ccl2", "Gata2", "Gata3", "Sox2", "Oct4", "Nanog", "Dab2"))
FeaturePlot(matrix.su, features = c( "Klf17", "Tfap2c", "Slc28a3", "Adap2", "Igfbp3", "Bamb1", "Havcr1"))
# Cga

#### DECIDE ON PARAMETERS:::
#l2fc_cutoff <- 0.6
l2fc_cutoff <- 0.25
database <- "KEGG_2019_Mouse"
# WikiPathways_2019_Mouse, GO_Biological_Process_2018, GO_Cellular_Component_2018, GO_Molecular_Function_2018,  BioCarta_2016, KEGG_2019_Mouse

GO_matrix <- data.frame(c("xx", "xxx"))
colnames(GO_matrix) <- paste0(database,".Term")  

for (z in seq_along(result_list)){
  res <- result_list[[z]]
  res_name <- names(result_list)[[z]]
  for (w in 1:length(res)){
    GO_matrix <- save_GO_res(w, res, db = database ) 
  }
}

dim(GO_matrix) # 
GO_matrix[,1] <- gsub(",",";", GO_matrix[,1])

#WP_matrix_l2fc0.6  <- GO_matrix
#GOBP_matrix_l2fc0.6  <- GO_matrix
#GOCC_matrix_l2fc0.6  <- GO_matrix
#GOMF_matrix_l2fc0.6  <- GO_matrix
#BioCarta_matrix_l2fc0.6  <- GO_matrix
#Kegg_matrix_l2fc0.6 <- GO_matrix

#WP_matrix_l2fc0.25 <- GO_matrix
#GOBP_matrix_l2fc0.25 <- GO_matrix
#GOCC_matrix_l2fc0.25 <- GO_matrix
#GOMF_matrix_l2fc0.25 <- GO_matrix
#BioCarta_matrix_l2fc0.25 <- GO_matrix
#Kegg_matrix_l2fc0.25 <- GO_matrix

GO_res_list <- list(GOBP_matrix_l2fc0.25=GOBP_matrix_l2fc0.25,
                    GOCC_matrix_l2fc0.25=GOCC_matrix_l2fc0.25,
                    GOMF_matrix_l2fc0.25=GOMF_matrix_l2fc0.25,
                    BioCarta_matrix_l2fc0.25=BioCarta_matrix_l2fc0.25,
                    WP_matrix_l2fc0.25=WP_matrix_l2fc0.25,
                    Kegg_matrix_l2fc0.25=Kegg_matrix_l2fc0.25)
#GO_res_list <- list(GOBP_matrix_l2fc0.6=GOBP_matrix_l2fc0.6,
#                    GOCC_matrix_l2fc0.6=GOCC_matrix_l2fc0.6,
#                    GOMF_matrix_l2fc0.6=GOMF_matrix_l2fc0.6,
#                    BioCarta_matrix_l2fc0.6=BioCarta_matrix_l2fc0.6,
#                    WP_matrix_l2fc0.6=WP_matrix_l2fc0.6,
#                    Kegg_matrix_l2fc0.6=Kegg_matrix_l2fc0.6)

#saveRDS(GO_res_list, "NewAlign51_SCT_HARMONY_res1.0_GO_res_list_l2fc0.25.Rds")
names(GO_res_list)

for(i in seq_along(GO_res_list)){
  GO_res_list[[i]] <- GO_res_list[[i]][-c(1,2),]
}
head(GO_res_list[[1]],2)

for(i in seq_along(GO_res_list)){
  write.csv(GO_res_list[[i]], paste( "NewAlign51_SCT_HARMONY_res_1.0_", names(GO_res_list)[[i]], ".csv", sep = "_"),  quote = FALSE)
}




#names(GO_res_list)[6] <-"Kegg_matrix_l2fc0.6" 
#names(GO_res_list)[12] <-"Kegg_matrix_l2fc0.25"

GO_matrix_for_plotting <- GO_res_list[[4]]
GO_res_name <- names(GO_res_list)[[4]]

for(i in seq_along(GO_res_list)){
#for(i in c(7,8)){
  GO_matrix_for_plotting <- GO_res_list[[i]]
  GO_res_name <- names(GO_res_list)[[i]]
  l2fc_cutoff <- as.numeric(as.character(gsub(".*l2fc", "", GO_res_name)))
  GO_matrix_for_plotting <- GO_matrix_for_plotting[GO_matrix_for_plotting[,1] != c("xx", "xxx"),]
  GO_matrix_molten <- reshape2::melt(GO_matrix_for_plotting)
  GO_matrix_molten <- na.omit(GO_matrix_molten)
  GO_matrix_molten <- as.data.frame(table(GO_matrix_molten[,1]))
  GO_matrix_molten <- GO_matrix_molten[order(GO_matrix_molten$Freq, decreasing = T),]
  
  GO_matrix <- unique(GO_matrix_for_plotting)
  GO_matrix2 <- GO_matrix
  GO_matrix2[1:5,1:5]
  rownames(GO_matrix2) <- GO_matrix2[,1]
  GO_matrix2 <- GO_matrix2[,-1]
  
  GO_matrix2_clust1 <- colnames(GO_matrix2)
  GO_matrix2_clust1 <- gsub(".vs.*", "", GO_matrix2_clust1)
  GO_matrix2_clust2 <- colnames(GO_matrix2)
  GO_matrix2_clust2 <- gsub(".*.vs.", "", GO_matrix2_clust2)
  GO_matrix2_clust2 <- gsub("X", "", GO_matrix2_clust2)
  
  rownames(GO_matrix2) <- gsub( " \\(GO:.*" ,  "" , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "regulation" ,  "reg." , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "dependent" ,  "dep." , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "positive" ,  "+" , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "negative" ,  "-" , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "pathway" ,  "path." , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "chemical" ,  "chem." , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "polymerase" ,  "pol" , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "cotranslational" ,  "cotransl." , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "biosynthetic" ,  "biosynth." , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "mitochondrial" ,  "MT" , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( "calcium" ,  "Ca" , rownames(GO_matrix2))
  rownames(GO_matrix2) <- gsub( " Homo sapiens h.*" ,  "" , rownames(GO_matrix2))
  
  dim(GO_matrix2)
  f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
  ht1 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", width = unit(ncol(GO_matrix2), "cm"),) # width = unit(140, "cm"),
  ht1
  
  pdf(paste( Project, "ComplexHeatmap", "Fig__ALL_Pairwise_", GO_res_name, reso, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix2)/2+5), height=nrow(GO_matrix2)/2)
  par(bg=NA)
  draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
}





