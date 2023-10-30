#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()
options(bitmapType='cairo')


library(cowplot)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(clustree)
library(harmony)
library(viridis)
library(colorspace)
library(optrees)
library(deldir)
#library(alphahull)

options(future.globals.maxSize = 4000 * 1024^2)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
Project <- "CTR_dscj1_HARMONY_orig_400g"




dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_clust24_25rm.Rds")
nrow(dropseq.integrated@meta.data[dropseq.integrated@meta.data$orig.ident2 == "Batch0004_9617_N721",])
nrow(dropseq.integrated@meta.data[dropseq.integrated@meta.data$orig.ident2 == "Batch0004_9617_N729",])


#saveRDS(matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/matrix.su_400g_SCT.Rds")



message("+-------------------------------------------------------------------------------")
message("+               load in dataset for harmony integration                         ")
message("+-------------------------------------------------------------------------------")


matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")
matrix.su <- readRDS("matrix.su_400g_SCT.Rds") # tsting if 400g or 500g cutoff works better???
matrix.su <- readRDS("matrix.su_500g_SCT.Rds") # tsting if 400g or 500g cutoff works better???


nrow(matrix.su@meta.data)             # 229793 
min(matrix.su@meta.data$nFeature_RNA) # 300
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident == "Batch0004_9617_N721",])
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident == "Batch0004_9617_N729",])
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N721")
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N729")

head(matrix.su@meta.data,2)
tmp_df <- unique(matrix.su@meta.data[,c("Experiment", "orig.ident2", "Batch")])
table(tmp_df[,c(1,3)])

matrix.su <- subset(matrix.su,  nFeature_RNA > 400) # 169094 cells
min(matrix.su@meta.data$nFeature_RNA)


message("+-------------------------------------------------------------------------------")
message("+                                HARMONY orig                               ")
message("+-------------------------------------------------------------------------------")

harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", reduction.save = "harmony") 
gc()
harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony", dims = 1:30) 
gc()
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.6)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.8)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1.5)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 2)

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 25, repulsion.strength = 2, spread = 2L )

#saveRDS(harmony_matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")

#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_clust24_25rm.Rds")




message("+-------------------------------------------------------------------------------")
message("+                                HARMONY lambda0.6                              ")
message("+-------------------------------------------------------------------------------")

harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", lambda = 0.6, reduction.save = "harmony_lambda0.6") # more diverse clusters

harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony_lambda0.6", dims = 1:30) 
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1.5)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.6)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.8)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0)

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony_lambda0.6", dims = 1:30, min.dist = 0.01)

#saveRDS(harmony_matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_lambda0.6.Rds")



message("+-------------------------------------------------------------------------------")
message("+                                HARMONY sigma0.5                               ")
message("+-------------------------------------------------------------------------------")

harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", sigma = 0.5, reduction.save = "harmony_sigma0.5") # more diverse clusters

harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony_sigma0.5", dims = 1:30) 
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1.5)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.6)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.8)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0)

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony_sigma0.5", dims = 1:30, min.dist = 0.01)

#saveRDS(harmony_matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_sigma0.5.Rds")


message("+-------------------------------------------------------------------------------")
message("+                                HARMONY STD                               ")
message("+-------------------------------------------------------------------------------")

harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", reduction.save = "harmony") # more diverse clusters

harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony", dims = 1:30) 
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1.5)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.6)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.8)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0)

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01)

#saveRDS(harmony_matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_STD.Rds")



message("+-------------------------------------------------------------------------------")
message("+                                HARMONY theta4                                 ")
message("+-------------------------------------------------------------------------------")

harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", theta = 4, reduction.save = "harmony_theta4") # more diverse clusters

harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony_theta4", dims = 1:30) 
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1.5)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.6)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.8)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0)

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony_theta4", dims = 1:30, min.dist = 0.01)

#saveRDS(harmony_matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_theta4.Rds")

dropseq.integrated <- RunUMAP(dropseq.integrated, reduction = "harmony_theta4", dims = 1:30, min.dist = 0.01, n.neighbors = 25, repulsion.strength = 2, spread = 2L )

#saveRDS(dropseq.integrated, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_theta4_repuls2_spread2.Rds")




############### remove for the time being....
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





message("+-------------------------------------------------------------------------------")
message("+                    load in  HARMONY object                                    ")
message("+-------------------------------------------------------------------------------")

#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_STD.Rds")
#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_theta4.Rds")
#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_sigma0.5.Rds")
#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_harmony_lambda0.6.Rds")
dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729.Rds")
dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_400g.Rds")

Project <- "CTR_dscj1_HARMONY_orig_400g"
gc()

head(dropseq.integrated@meta.data,2)



message("+-------------------------------------------------------------------------------")
message("+              Explore and QC integrated dataset                                ")
message("+-------------------------------------------------------------------------------")
  
plots <- DimPlot(dropseq.integrated, group.by = c("Batch", "Experiment"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))


plots2 <- DimPlot(dropseq.integrated, group.by = "Batch", split.by = "Experiment", ncol = 3)



message("+-------------------------------------------------------------------------------")
message("+                       CLUSTREE                                                ")
message("+-------------------------------------------------------------------------------")

head(dropseq.integrated@meta.data,2)
#dropseq.integrated@meta.data <- dropseq.integrated@meta.data[,-c(18,19)]

tree_to_plot <- clustree::clustree(dropseq.integrated, prefix = "SCT_snn_res.") # "RNA_snn_res."

pdf(paste(Project, "harmony_400g___matrix.su_clustree_0_to_2.pdf", sep=""), width=12,height=8)
par(bg=NA)
tree_to_plot
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ ONLY ONCE DECIDE ON RESOLUTION: remove clusters smaller than 100 cells        ")
message("+-------------------------------------------------------------------------------")

# clusters 24 and 25 have 34 and 31 cells respectively.
head(dropseq.integrated@meta.data,2)
clusts_below100 <- as.data.frame(table(dropseq.integrated@meta.data[,"SCT_snn_res.1"]))
keep_clusters_over_100 <- clusts_below100[clusts_below100$Freq >= 100,]$Var1

dropseq.integrated <- subset(dropseq.integrated, SCT_snn_res.1 %in% keep_clusters_over_100)


message("+-------------------------------------------------------------------------------")
message("+        dot plot of TSC markers                                 ")
message("+-------------------------------------------------------------------------------")

head(dropseq.integrated@meta.data,2)
nrow(dropseq.integrated@meta.data) # after removing clusters < 100 cells :  216547

dropseq.integrated@meta.data$Cluster <- as.character(dropseq.integrated@meta.data$SCT_snn_res.1)
table(dropseq.integrated@meta.data$Cluster)

dropseq.integrated@meta.data$Cluster <- factor(dropseq.integrated@meta.data$Cluster , levels = c("12", "6", "18", "1", "19", "10", "23", "16", "2", "14", "11", "5", "15", "0", "9", "17", "8", "13", "4", "20", "3", "21", "7","22"))

Idents(dropseq.integrated) <- (dropseq.integrated@meta.data$Cluster)


sum(dropseq.integrated@meta.data$nFeature_RNA > 500)

#DefaultAssay(dropseq.integrated)
DotPlot(dropseq.integrated, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18") , cols = c("green", "blue"), assay = "RNA", dot.scale = 10, scale = TRUE)

DotPlot(dropseq.integrated, features = c("Tfeb", "Ly6e","Ovol2","Gcm1","Dlx3","Cebpa", "Syna") , cols = c("green", "blue"), assay = "RNA", dot.scale = 10, scale = TRUE)
DotPlot(dropseq.integrated500, features = c("Tfeb", "Ly6e","Ovol2","Gcm1","Dlx3","Cebpa", "Syna") , cols = c("green", "blue"), assay = "RNA", dot.scale = 10, scale = TRUE)


message("+-------------------------------------------------------------------------------")
message("+                        Plotting UMAP - custom                                 ")
message("+-------------------------------------------------------------------------------")

sampleTable <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/sampleTable_MERGED_53.csv")

#dataset_name <- "Harmony__STD" 
dataset_name <- "Harmony__orig" 
#dataset_name <- "Harmony__theta4_repuls2spread2" 
#dataset_name <- "Harmony__sigma0.5" 
#dataset_name <- "Harmony__lambda0.6" 

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
  Cluster_Anno[i] <-  paste0(  "annotate(\'text\', label = paste0(\'", Clusters[i], "\'), x =", coord_x[i], ", y = ", coord_y[i], ", color = \'black\', size=8)" )
}

Cluster_Anno <- unlist(Cluster_Anno)

library(stringr)
Cluster_Anno <- str_c(Cluster_Anno,  collapse =  " + ")


centers <- data.frame(Clusters=Clusters,coord_x=unlist(coord_x), coord_y=unlist(coord_y))
centers <- centers[order(centers$Clusters),]

dat <- matrix.umap

message("----------------------go to script ASSIGN_COLOURS_UMAP.R !!! -----------------------")

gc()

message("----------------------umap.Cluster-----------------------")

unique(matrix.umap$Cluster)
matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22","23")) 

Col_cluster <- c(best.colordf$color)
names(Col_cluster) <- levels(matrix.umap$Cluster)

# for harmony orig after removing lcust 24 and 25:::
Col_cluster <- c("#cc99ff", "#0868ac" ,"#8c96c6" ,  "#FF9900FF", "#810f7c",
                 "#004529",   "#4d004b", "#333399"  , "#00B3FFFF" ,"#BF5B17", 
                 "#FF00E6FF" ,"#7F00FFFF", "#FFFF99" ,  "#f768a1"  , "#00FFB2FF",
                 "#78c679"   ,  "#d9f0a3" ,  "#bfd3e6" , "#88419d" ,  "#666666"  , 
                 "#d4b9da"  , "#fa9fb5"  , "#FF0000FF" ,"#fde0dd")
names(Col_cluster) <- levels(matrix.umap$Cluster)

# for harmony std
Col_cluster<-c("#00B3FFFF","#8c96c6","#78c679","#d4b9da","#d9f0a3","#fde0dd","#BF5B17","#fa9fb5","#810f7c","#004529","#4d004b","#f768a1","#4eb3d3","#FF0000FF","#081d58","#FF00E6FF","#FFFF99","#FF9900FF","#FFE500FF","#88419d","#7F00FFFF","#00FFB2FF","#238443","#0868ac","#bfd3e6")
names(Col_cluster) <- levels(matrix.umap$Cluster)

# for harmony theta4
Col_cluster<-c("#00B3FFFF", "#4d004b" ,  "#FFE500FF" ,"#FFFF99" ,  "#FF0000FF", "#FF00E6FF","#BF5B17"  ,"#810f7c" ,    "#004529" ,  "#f768a1" ,  "#081d58" ,  "#FF9900FF", "#d9f0a3" ,  "#88419d" ,  "#238443" ,  "#fa9fb5" ,  "#78c679" ,  "#7F00FFFF" ,"#00FFB2FF" ,"#d4b9da",  "#fde0dd" ,  "#0868ac",   "#8c96c6" ,  "#4eb3d3" ,  "#bfd3e6")
names(Col_cluster) <- levels(matrix.umap$Cluster)


# for harmony lambda0.6
Col_cluster <- c("#666666" ,  "#d4b9da" ,  "#FF9900FF" ,"#00B3FFFF", "#BF5B17" ,  "#fa9fb5"  , "#FF0000FF" ,"#004529" ,  "#4d004b" ,  "#00FFB2FF", "#081d58" ,  "#f768a1",   "#88419d" ,  "#810f7c",   "#7F00FFFF" ,"#8c96c6"  , "#FF00E6FF", "#0868ac" ,  "#FFE500FF", "#bfd3e6" ,"#fde0dd" ,  "#4eb3d3" ,  "#238443" ,  "#FFFF99" ,  "#d9f0a3" ,  "#78c679")
names(Col_cluster) <- levels(matrix.umap$Cluster)

# for harmony sigma0.5
Col_cluster <- c("#4eb3d3" ,  "#f768a1" ,  "#FF0000FF", "#004529" ,  "#78c679" ,  "#fde0dd" ,  "#BF5B17" ,  "#238443",   "#FFFF99"  , "#7F00FFFF","#FF9900FF" ,"#FFE500FF", "#d9f0a3",   "#00FFB2FF", "#bfd3e6"  , "#d4b9da"  , "#00B3FFFF" ,"#fa9fb5",   "#666666" ,  "#8c96c6" ,"#081d58" ,  "#0868ac",   "#4d004b" ,  "#810f7c" ,  "#FF00E6FF", "#333399"  , "#88419d")
names(Col_cluster) <- levels(matrix.umap$Cluster)

# for harmony theta4 - repuls2 spread2:
#"#fa9fb5"   "#d9f0a3"   "#4eb3d3"   "#238443"   "#fde0dd"   "#00FFB2FF"
#[7] "#78c679"   "#004529"   "#d4b9da"   "#810f7c"   "#FF0000FF" "#FFE500FF"
#[13] "#00B3FFFF" "#f768a1"   "#cc99ff"   "#FF9900FF" "#0868ac"   "#bfd3e6"  
#[19] "#8c96c6"   "#4d004b"   "#081d58"   "#FF00E6FF" "#333399"   "#88419d"  
#[25] "#7F00FFFF"

# for harmony orig:::
Col_cluster<-c( "#FF00E6FF", "#FF9900FF", "#FF0000FF", "#238443" , "#fa9fb5" , 
                "#00B3FFFF",   "#88419d" ,  "#810f7c" ,  "#081d58" , "#FFE500FF",
                "#FFFF99"  ,  "#f768a1" ,  "#004529" ,  "#8c96c6" ,  "#4d004b" ,
                "#fde0dd" ,  "#7F00FFFF", "#d9f0a3" ,  "#0868ac"  ,"#d4b9da"  ,
                "#4eb3d3" ,  "#bfd3e6",   "#78c679" ,   "#00FFB2FF")
names(Col_cluster) <- levels(matrix.umap$Cluster)


#Cluster_Anno

umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  #ggtitle(paste0(dataset_name, " QC :::  ", reso, " Cluster" )) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=1, size=6))) + 
  theme_classic() + 
  scale_colour_manual("Cluster", values = Col_cluster) +
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) #+ theme(legend.position = "none")

#umap.Cluster_w_legend <- umap.Cluster
#umap.Cluster_no_legend <- umap.Cluster

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Cluster_w_legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Cluster_no_legend
dev.off()

#saveRDS(umap.Cluster_w_legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.Cluster_w_legend.Rds")
#saveRDS(umap.Cluster_no_legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.Cluster_no_legend.Rds")

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Cluster_w_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Cluster_no_legend
dev.off()


message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  #ggtitle(paste0( dataset_name ,"QC ::: ", reso, " Cell Cycle" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() #+ theme(legend.position = "none") 

#umap.cc_w_Legend <- umap.cc
#umap.cc_no_Legend <- umap.cc


png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.cc_w_Legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.cc_no_Legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.cc_w_Legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.cc_no_Legend
dev.off()


#saveRDS(umap.cc_w_Legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.cc_w_Legend.Rds")
#saveRDS(umap.cc_no_Legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.cc_no_Legend.Rds")


message("----------------------umap.Experiment-----------------------")

umap.Experiment    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.2, size=0.3) +
  #ggtitle(paste0("Resolution ", reso, " :: Experiment " )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Experiment, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="black", "I.1"="yellow", "R.1"="yellow2", 
                                     "I.4"="orangered", "R.4"="red", "I.24"="lightblue2", "R.24"="lightblue4",
                                     "I.36"="blue", "R.36"="darkblue","I.48"="purple", "R.48"="orchid"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  theme_classic() + coord_fixed() +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=6), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) #+ theme(legend.position = "none") 

#umap.Expt_w_Legend <- umap.Experiment
#umap.Expt_no_legend <- umap.Experiment


png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Expt_w_legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Expt_no_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Expt_w_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Expt_no_legend
dev.off()


#saveRDS(umap.cc_w_Legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.cc_w_Legend.Rds")
#saveRDS(umap.cc_no_Legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.cc_no_Legend.Rds")


Cluster_Anno

message("----------------------umap.t.0-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_timepoints")
table(matrix.umap$Mixer)

umap.t0    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.3, size=0.3) +
  #ggtitle(paste0("Resolution ", reso, " T0 in Batches" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "0.0" = "black")) +
  coord_fixed() +
  xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  #guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + theme(legend.position = "none")





message("----------------------umap.t.1-----------------------")

umap.Experiment1    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="blue", "R.1"="red", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + theme(legend.position = "none")





message("----------------------umap.t.4-----------------------")

umap.Experiment4    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="blue", "R.4"="red", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + theme(legend.position = "none")






message("----------------------umap.t.24-----------------------")

umap.Experiment24    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="blue", "R.24"="red",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c('0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))  + 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + theme(legend.position = "none")





message("----------------------umap.t.36-----------------------")

umap.Experiment36    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="blue", "R.36"="red","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + theme(legend.position = "none")



message("----------------------umap.t.48-----------------------")

umap.Experiment48    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="blue", "R.48"="red"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) + theme(legend.position = "none")




png(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "umap__Expt_grid_eachTimePoint", reso, ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
plot_grid(umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 2, align = "hv")
dev.off()


#####
list_umaps_for_plotting_noLegend <- list(umap.t0=umap.t0, 
                                         umap.Experiment1=umap.Experiment1, 
                                         umap.Experiment4=umap.Experiment4, 
                                         umap.Experiment24=umap.Experiment24, 
                                         umap.Experiment36=umap.Experiment36, 
                                         umap.Experiment48=umap.Experiment48)
#saveRDS(list_umaps_for_plotting_noLegend,"list_umaps_for_plotting_noLegend.Rds")
list_umaps_for_plotting_2 <- list(umap.cc_no_Legend=umap.cc_no_Legend, 
                                  umap.cc_w_Legend =umap.cc_w_Legend, 
                                  umap.Cluster_no_legend  =umap.Cluster_no_legend, 
                                  umap.Cluster_w_legend=umap.Cluster_w_legend, 
                                  umap.Expt_no_legend =umap.Expt_no_legend, 
                                  umap.Expt_w_Legend=umap.Expt_w_Legend )
#saveRDS(list_umaps_for_plotting_2,"list_umaps_for_plotting2_cc_clust_expt.Rds")


list_umaps_for_plotting_noLegend <- readRDS("list_umaps_for_plotting_noLegend.Rds")
list_umaps_for_plotting_2        <- readRDS("list_umaps_for_plotting2_cc_clust_expt.Rds")
plt_pseudotime <- readRDS("MONOCLE_SCT/harmony-orig_monocle2_plt_pseudotime_loopT.Rds")

names(list_umaps_for_plotting_noLegend)
names(list_umaps_for_plotting_2)

pdf(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "9x_GRID_PAPER_FIGURE", reso, ".pdf", sep="_"), width=20,height=18)
#png(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "9x_GRID_PAPER_FIGURE", reso, ".png", sep="_"), width=2500,height=2000, type = "cairo")
par(bg=NA)
plot_grid( list_umaps_for_plotting_2[["umap.Cluster_no_legend"]], 
           list_umaps_for_plotting_2[["umap.cc_no_Legend"]] + theme(legend.position = "none")+ xlab("") + ylab("") + coord_fixed() ,
           plt_pseudotime + theme(legend.position = "none") + xlab("") + ylab("")+ coord_fixed() , 
           list_umaps_for_plotting_noLegend[["umap.t0"]], 
           list_umaps_for_plotting_noLegend[["umap.Experiment1"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment4"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment24"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment36"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment48"]], ncol = 3, align = "vh")
dev.off()


#need to calculate pseudotime!!! 
png(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "GRID_PAPER_FIGURE", reso, ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
plot_grid(umap.Cluster_no_legend, umap.cc_no_Legend, umap.pseudotime, list_umaps_for_plotting_noLegend["umap.t0"], list_umaps_for_plotting_noLegend["umap.Experiment1"],list_umaps_for_plotting_noLegend["umap.Experiment4"],list_umaps_for_plotting_noLegend["umap.Experiment24"],list_umaps_for_plotting_noLegend["umap.Experiment36"],list_umaps_for_plotting_noLegend["umap.Experiment48"], ncol = 3, nrow = 3, align = "hv")
dev.off()




#####pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "GRID_PAPER_FIGURE","res_", reso,".pdf"), width=20, height=10)
par(bg=NA)
plot_grid(umap.Cluster_no_legend, umap.cc_no_Legend, umap.pseudotime, umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 3, align = "hv")
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


unique(matrix.umap[matrix.umap$Experiment == "R.36",]$BatchSample)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36", "R.36", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R4_reps_N712", "0B0_Dups_36R4_reps_N712", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R2_reps_N707", "0B0_Dups_36R2_reps_N707", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R1_reps_N706", "0B0_Dups_36R1_reps_N706", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R3_reps_N710", "0B0_Dups_36R3_reps_N710", matrix.umap$QC)

umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC R.36 samples" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "0B0_Dups_36R4_reps_N712" = "green", "0B0_Dups_36R2_reps_N707"="purple", "0B0_Dups_36R1_reps_N706"="blue", "0B0_Dups_36R3_reps_N710"="red")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc



unique(matrix.umap[matrix.umap$Experiment == "R.4",]$BatchSample)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4", "R.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "ABC_Dups_5A_reps_N726", "ABC_Dups_5A_reps_N726", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "AB0_Dups_6C_reps_N720", "AB0_Dups_6C_reps_N720", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "ABC_Dups_6B_reps_N715", "ABC_Dups_6B_reps_N715", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "A0C_Dups_6D_reps_N720", "A0C_Dups_6D_reps_N720", matrix.umap$QC)

umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC R.4 samples" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_Dups_5A_reps_N726" = "green", "AB0_Dups_6C_reps_N720"="purple", "ABC_Dups_6B_reps_N715"="blue", "A0C_Dups_6D_reps_N720"="red")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc



unique(matrix.umap[matrix.umap$Experiment == "I.24",]$BatchSample)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24", "I.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "ABC_Dups_8D_reps_N722", "ABC_Dups_8D_reps_N722", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "ABC_Dups_7A_reps_N728", "ABC_Dups_7A_reps_N728", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "AB0_Dups_8C_reps_N706", "AB0_Dups_8C_reps_N706", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "ABC_Dups_8B_reps_N710", "ABC_Dups_8B_reps_N710", matrix.umap$QC)

umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC I.24 samples" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_Dups_8D_reps_N722" = "green", "ABC_Dups_7A_reps_N728"="purple", "AB0_Dups_8C_reps_N706"="blue", "ABC_Dups_8B_reps_N710"="red")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc










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
message("+                        FIND MARKERS for res 1.                                ")
message("+-------------------------------------------------------------------------------")

Project <- "CTR_dscj1_HARMONY_markers"

markerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_together_MARKERS/"
setwd(markerDir)

# set up variables:
l2fc_cutoff <- 0.25
reso <- "res.1"
last_cluster <- length(unique(dropseq.integrated@meta.data[, paste0("SCT_snn_", reso)]))-1

head(dropseq.integrated@meta.data,2)
Idents(dropseq.integrated) <-  dropseq.integrated@meta.data[,paste0("SCT_snn_", reso)]

cluster_numbers <- c(0:(  length(unique(dropseq.integrated@meta.data[, paste0("SCT_snn_", reso)]))-1 ))
marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = length(cluster_numbers), nrow = length(cluster_numbers) ))

marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = last_cluster+1, nrow = last_cluster+1 ))
rownames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20","c21","c22", "c23") 
colnames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20","c21","c22", "c23") 

table(dropseq.integrated@meta.data[, paste0("SCT_snn_", reso)])




library(purrr)

findMarkers_for_cluster_pair <- function(i){
  if (j == i){
    print("same cluster")
  } else {
    tmp_markers <- FindMarkers(matrix.su_sep, ident.1 = j, ident.2= i, logfc.threshold = l2fc_cutoff, verbose = TRUE)
    print(length(unique(row.names(tmp_markers))))
    return(tmp_markers)
  }
}
safely_findMarkers = safely( findMarkers_for_cluster_pair )


for (n in 0:last_cluster){
  cluster_numbers <- c(n:last_cluster)
  j <- as.character(n)
  result = map(cluster_numbers, safely_findMarkers)
  saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
  result_list[[n]] <- result
}
  

cluster_numbers <- c(0:last_cluster)
j <- "0"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_00 <- result

cluster_numbers <- c(1:last_cluster)
j <- "1"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_01 <- result

cluster_numbers <- c(2:last_cluster)
j <- "2"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_02 <- result

cluster_numbers <- c(3:last_cluster)
j <- "3"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name,reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_03 <- result

cluster_numbers <- c(4:last_cluster)
j <- "4"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_04 <- result


cluster_numbers <- c(5:last_cluster)
j <- "5"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_05 <- result

cluster_numbers <- c(6:last_cluster)
j <- "6"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_06 <- result

cluster_numbers <- c(7:last_cluster)
j <- "7"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_07 <- result

cluster_numbers <- c(8:last_cluster)
j <- "8"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_08 <- result

cluster_numbers <- c(9:last_cluster)
j <- "9"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso, j, "_result_", "_l2fc" , l2fc_cutoff ,  ".rds"))
result_09 <- result

cluster_numbers <- c(10:last_cluster)
j <- "10"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_10 <- result

cluster_numbers <- c(11:last_cluster)
j <- "11"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_11 <- result

cluster_numbers <- c(12:last_cluster)
j <- "12"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_12 <- result

cluster_numbers <- c(13:last_cluster)
j <- "13"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_13 <- result

cluster_numbers <- c(14:last_cluster)
j <- "14"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_14 <- result

cluster_numbers <- c(15:last_cluster)
j <- "15"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_15 <- result

cluster_numbers <- c(16:last_cluster)
j <- "16"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_16 <- result

cluster_numbers <- c(17:last_cluster)
j <- "17"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_17 <- result

cluster_numbers <- c(18:last_cluster)
j <- "18"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_18 <- result

cluster_numbers <- c(19:last_cluster)
j <- "19"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_19 <- result

cluster_numbers <- c(20:last_cluster)
j <- "20"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_20 <- result

cluster_numbers <- c(21:last_cluster)
j <- "21"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_21 <- result

cluster_numbers <- c(22:last_cluster)
j <- "22"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_22 <- result

cluster_numbers <- c(23:last_cluster)
j <- "23"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_23 <- result



result_list <- list(result_00  = result_00, 
                    result_01 = result_01, 
                    result_02 = result_02, 
                    result_03 = result_03, 
                    result_04 = result_04, 
                    result_05 = result_05,
                    result_06 = result_06,
                    result_07 = result_07,
                    result_08 = result_08,
                    result_09 = result_09,
                    result_10 = result_10,
                    result_11 = result_11,
                    result_12 = result_12,
                    result_13 = result_13,
                    result_14 = result_14,
                    result_15 = result_15,
                    result_16 = result_16,
                    result_17 = result_17,
                    result_18 = result_18,
                    result_19 = result_19,
                    result_20 = result_20,
                    result_21 = result_21,
                    result_22 = result_22,
                    result_23 = result_23)
names(result_list)










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





