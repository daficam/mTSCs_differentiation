#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()
options(bitmapType='cairo')


library(cowplot)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(sceasy)

options(future.globals.maxSize = 4000 * 1024^2) # limit 4Gb

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I"
setwd(baseDir)
Project <- "CTR_dscj1_harmony_R_I_sep"

sampleTable <- read.csv("sampleTable_MERGED_53.csv")


Project <- "CTR_dscj1_harmony_R_500g"
matrix.su_sep <- readRDS("HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_R.Rds")
#dataset_name <- "Remove_29_500g"
matrix.umap <- readRDS("Harmony_Remove_29_500g_matrix.umap_.Rds")

Project <- "CTR_dscj1_harmony_I_500g"
matrix.su_sep <- readRDS("HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_I.Rds")
dataset_name <- "Inhibit_32_500g"
matrix.umap <- readRDS("Harmony_Inhibit_32_500g_matrix.umap_.Rds")


tmp_tbl <- as.data.frame(table(matrix.su_sep@meta.data[,c(reso, "Experiment")]))
#write.csv(tmp_tbl_I, "Harmony_I_500g_TABLE_no_of_cells_per_Cluster_Experiment.csv")

tmp_tbl_I <- reshape2::dcast(tmp_tbl_I, formula = SCT_snn_res.0.8 ~ Experiment)
tmp_tbl_R <- reshape2::dcast(tmp_tbl_R, formula = SCT_snn_res.0.8 ~ Experiment)
tmp_tbl_orig <- reshape2::dcast(tmp_tbl_orig, formula = SCT_snn_res.0.8 ~ Experiment)
#write.csv(tmp_tbl_I, "Harmony_I_500g_TABLE_no_of_cells_per_Cluster_Experiment.csv")
#write.csv(tmp_tbl_R, "Harmony_R_500g_TABLE_no_of_cells_per_Cluster_Experiment.csv")
#write.csv(tmp_tbl_orig, "Harmony_orig_500g_TABLE_no_of_cells_per_Cluster_Experiment.csv")


#   Seurat to AnnData
#sceasy::convertFormat(matrix.su_sep, from= "seurat", to= "anndata", outFile = "matrix.su_sep_Remove_29_500g.h5ad")
#sceasy::convertFormat(matrix.su_sep, from= "seurat", to= "anndata", outFile = "matrix.su_sep_Inhibit_32_500g.h5ad")


# FeaturePlot: Phf8, Gcm1, E2f8, ...
pdf(paste( "FeaturePlot", dataset, "_green_Phf8", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
FeaturePlot(matrix.su_sep,  features =c("Phf8","Gcm1","E2f8"), slot = "data", ncol = 3, pt.size =0.3 , max.cutoff  =  1.5)
dev.off()



message("--------------------------------------------------------------------------------")
message("+           RUN HARMONY ON SEPARATE R AND I OBJECTS                             ")
message("+-------------------------------------------------------------------------------")

#matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53.Rds")
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N721",])
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N729",])


matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N729")
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N721")

matrix.su_I <- subset(matrix.su, Treatment != "R")
matrix.su_R <- subset(matrix.su, Treatment != "R")

matrix.su_I <- PercentageFeatureSet(matrix.su_I, pattern = "^mt-", col.name = "percent.mt")
matrix.su_R <- PercentageFeatureSet(matrix.su_R, pattern = "^mt-", col.name = "percent.mt")
max(matrix.su_R@meta.data$percent.mt) # 73.59477%
matrix.su_R <- subset(matrix.su_R, subset = percent.mt < 5)
matrix.su_I <- subset(matrix.su_I, subset = percent.mt < 5)

matrix.su_I <- SCTransform(matrix.su_I, return.only.var.genes = FALSE, vars.to.regress = c("percent.mt"), verbose = T)
matrix.su_R <- SCTransform(matrix.su_R, return.only.var.genes = FALSE, vars.to.regress = c("percent.mt"), verbose = T)

matrix.su_I <- RunPCA(matrix.su_I, assay = "SCT", npcs = 50)
matrix.su_R <- RunPCA(matrix.su_R, assay = "SCT", npcs = 50)

matrix.su_I <- RunUMAP(matrix.su_I, reduction = "pca", assay = "SCT", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 
matrix.su_R <- RunUMAP(matrix.su_R, reduction = "pca", assay = "SCT", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 

length(unique(matrix.su_I@meta.data$orig.ident))
length(unique(matrix.su_R@meta.data$orig.ident))

head(matrix.su_I@meta.data,2)
matrix.su_I <- RunHarmony(matrix.su_I, group.by.vars = "Batch", reduction = "pca", assay.use="SCT")
matrix.su_R <- RunHarmony(matrix.su_R, group.by.vars = "Batch", reduction = "pca", assay.use="SCT")
gc()
matrix.su_I <- RunUMAP(matrix.su_I, reduction = "harmony", dims = 1:30, min.dist = 0.01)
matrix.su_R <- RunUMAP(matrix.su_R, reduction = "harmony", dims = 1:30, min.dist = 0.01)
gc()
matrix.su_I <- FindNeighbors(matrix.su_I, reduction = "harmony", dims = 1:30) 
matrix.su_R <- FindNeighbors(matrix.su_R, reduction = "harmony", dims = 1:30) 
gc()
matrix.su_I <- FindClusters(matrix.su_I, resolution = 0)
matrix.su_I <- FindClusters(matrix.su_I, resolution = 0.2)
matrix.su_I <- FindClusters(matrix.su_I, resolution = 0.6)
matrix.su_I <- FindClusters(matrix.su_I, resolution = 0.8)
matrix.su_I <- FindClusters(matrix.su_I, resolution = 1)
matrix.su_I <- FindClusters(matrix.su_I, resolution = 1.5)
matrix.su_I <- FindClusters(matrix.su_I, resolution = 2)
gc()
matrix.su_R <- FindClusters(matrix.su_R, resolution = 0)
matrix.su_R <- FindClusters(matrix.su_R, resolution = 0.2)
matrix.su_R <- FindClusters(matrix.su_R, resolution = 0.6)
matrix.su_R <- FindClusters(matrix.su_R, resolution = 0.8)
matrix.su_R <- FindClusters(matrix.su_R, resolution = 1)
matrix.su_R <- FindClusters(matrix.su_R, resolution = 1.5)
matrix.su_R <- FindClusters(matrix.su_R, resolution = 2)

gc()
matrix.su_R <- RunUMAP(matrix.su_R, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 25, repulsion.strength = 2, spread = 2L )
#saveRDS(matrix.su_R, "harmony_matrix.su_R__29__umap_n.n25_repuls2_spr2L.Rds")
matrix.su_I <- RunUMAP(matrix.su_I, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 25, repulsion.strength = 2, spread = 2L )
#saveRDS(matrix.su_I, "harmony_matrix.su_I__32__umap_n.n25_repuls2_spr2L.Rds")



message("--------------------------------------------------------------------------------")
message("+                    post HARMONY processing / qc                               ")
message("+-------------------------------------------------------------------------------")

setwd(ResDir)

#matrix.su_sep <- readRDS("harmony_matrix.su_R__29__umap_n.n25_repuls2_spr2L.Rds")
#dataset_name <- "Removal_29"

matrix.su_sep <- readRDS("harmony_matrix.su_I__32__umap_n.n25_repuls2_spr2L.Rds")
dataset_name <- "Inhibit_32"

gc()

head(matrix.su_sep@meta.data,2)
#matrix.su_sep <- FindNeighbors(matrix.su_sep, reduction = "harmony", dims = 1:30) 
#matrix.su_sep <- FindClusters(matrix.su_sep, resolution = 0.8)

message("--------------------------------------------------------------------------------")
message("+                                  cell cycle                                   ")
message("+-------------------------------------------------------------------------------")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- tolower(s.genes)
g2m.genes <- tolower(g2m.genes)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

s.genes <- firstup(s.genes)
g2m.genes <- firstup(g2m.genes)

matrix.su_sep <- CellCycleScoring(matrix.su_sep, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



message("+-------------------------------------------------------------------------------")
message("+                                 Clustree                                      ")
message("+-------------------------------------------------------------------------------")

tree_to_plot <- clustree::clustree(matrix.su_sep, prefix = "SCT_snn_res.") # "RNA_snn_res."

pdf(paste("matrix.su_REALIGNED_Harmony_separate_", dataset_name, "_300g___clustree.pdf", sep=""), width=12,height=8)
par(bg=NA)
tree_to_plot
dev.off()

colnames(matrix.su_sep@meta.data)
unique(matrix.su_sep@meta.data$SCT_snn_res.1) 
unique(matrix.su_sep@meta.data$SCT_snn_res.0.8)
unique(matrix.su_sep@meta.data$SCT_snn_res.0.6) 




message("+-------------------------------------------------------------------------------")
message("+        remove clusters smaller than 100 cells                                 ")
message("+-------------------------------------------------------------------------------")

# clusters 24 and 25 have 34 and 31 cells respectively.
head(matrix.su_sep@meta.data,2)
clusts_below100 <- as.data.frame(table(matrix.su_sep@meta.data[,"SCT_snn_res.0.8"]))
keep_clusters_over_100 <- clusts_below100[clusts_below100$Freq >= 100,]$Var1

matrix.su_sep <- subset(matrix.su_sep, SCT_snn_res.0.8 %in% keep_clusters_over_100)

#matrix.su_sep@meta.data$SCT_snn_res.0.8 <- factor(matrix.su_sep@meta.data$SCT_snn_res.0.8, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18" ))

table(matrix.su_sep@meta.data$SCT_snn_res.0.8)



message("+-------------------------------------------------------------------------------")
message("+      QC DOT PLOTS                                 ")
message("+-------------------------------------------------------------------------------")
# ORIG::: 300g:;:
#  R :::: 16 18 13 1  2  20 6  14 0  3  15 21 11 5  9  4  17 10 8  7  12 19
#  I :::: 9  6  1  17 14 8  11 0  7  3  4  12 13 10 18 15 2  5  16

Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.0.8

###### INHIBIT:

# from tempora for I 500g::: 
Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels= c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0"))

# from MONOCLE PSEUDOTIME for I 500g::: 
Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels= c("9",  "11" ,"16", "5" , "1" , "15", "6" , "12", "2" , "8" , "0" , "7" , "10", "14", "13", "4" , "3"))



###### REMOVAL:

# from tempora for R 500g::: 
Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels= c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8"))

# from MONOCLE PSEUDOTIME for R 500g::: 
Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels= c("14", "13", "2" , "3" , "7",  "0",  "15", "6" , "12", "11", "5" , "10", "9" , "4"  ,"8"  ,"1"))




# together??
#Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.1
#Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels= c("9",  "10", "7", "13", "1",  "17", "11", "15", "5",  "12",  "18",  "4",  "2", "6",  "14", "16" ,"8", "3","0")) # tempora !




# TSC sig
DotPlot(matrix.su_sep, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18") , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# EPC sig
DotPlot(matrix.su_sep, features = c("Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Gjb3", "Ascl2", "Tpbpa", "Prdm1", "Ets2", "Stra13", "Mdfi") , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #"Prl3b1", "Ctsq", 

# LP sig
DotPlot(matrix.su_sep, features = c("Tfeb", "Ly6e","Syna", "Gcm1", "Dlx3", "Cebpa", "Ovol2", "Junb", "Fosl1", "Arnt", "Hif1a", "Pparg") , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) #"Hif1b",,"Pparbp"


# ALL sig
DotPlot(matrix.su_sep, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt") , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(paste("DotPlot_", Project, dataset_name, "_", "umap.Cluster_no_legend_","res_", reso,".pdf", sep = "_"), width=9, height=4)
par(bg=NA)
DotPlot(matrix.su_sep, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




FeaturePlot(matrix.su_sep, features = c("Tfeb", "Ly6e","Syna", "Gcm1", "Dlx3", "Cebpa", "Ovol2", "Junb", "Fosl1", "Arnt", "Hif1a", "Pparg"))








message("+-------------------------------------------------------------------------------")
message("+                        Plotting UMAP - custom                                 ")
message("+-------------------------------------------------------------------------------")


head(matrix.su_sep@meta.data,2)
unique(matrix.su_sep@meta.data$SCT_snn_res.0.8)

matrix.umap            <- as.data.frame(Embeddings(object=matrix.su_sep, reduction="umap"))
matrix.umap$Sample_2   <- rownames(matrix.umap) 
matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
rownames(sampleTable)  <- sampleTable$Sample_2 
matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
matrix.umap$Project    <- sampleTable[matrix.umap$Sample_2, ]$Project
head(matrix.umap)

cell_cycle             <- matrix.su_sep@meta.data %>% dplyr::select((Phase))  # 
colnames(cell_cycle)   <- c("Phase")
matrix.umap$Phase      <- cell_cycle$Phase

reso                <- "SCT_snn_res.0.8"
clust               <- matrix.su_sep@meta.data %>% dplyr::select((SCT_snn_res.0.8))  # 

colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 

head(matrix.umap,2)
unique(matrix.umap$Batch)
unique(matrix.umap$Cluster)

unique(matrix.umap$Cluster)


# for *** res 0.8 ***
if(dataset_name == "Removal_29"){
  matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15" )) 
} else if(dataset_name == "Inhibit_32") {
  matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16" ))
} else { print("no dataset name chosen!!")}




table(matrix.umap$Experiment)
table(matrix.umap$Cluster)


message("----------------------Cluster annotation -----------------------")

#matrix.umap$UMAP_1 <- matrix.umap$UMAP_1 * -1
#matrix.umap$UMAP_2 <- matrix.umap$UMAP_2 * -1


Clusters <- as.numeric(levels(matrix.umap$Cluster))

coord_x <- list()
coord_y <- list()

for (i in 1:length(Clusters)){
  coord_x[i] <- mean(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$UMAP_1)  
  coord_y[i] <- mean(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$UMAP_2)  
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



message("----------------------umap.Cluster COLOURS -----------------------")

centers <- data.frame(Clusters=Clusters,coord_x=unlist(coord_x), coord_y=unlist(coord_y))
centers <- centers[order(centers$Clusters),]

dat <- matrix.umap

Col_cluster <- c(best.colordf$color)
names(Col_cluster) <- levels(matrix.umap$Cluster)


# for res 0.8 R only :::
Col_cluster<-c( "#d9f0a3" ,  "#bfd3e6" ,  "#FF9900FF" ,"#FFE500FF", "#f768a1" ,  "#4d004b" ,  "#00B3FFFF", "#4eb3d3" ,  "#78c679" ,  "#8c96c6" ,  "#fa9fb5" ,  "#00FFB2FF", "#0868ac"  , "#d4b9da",   "#fde0dd"  , "#238443" )
names(Col_cluster) <- levels(matrix.umap$Cluster)


# for res 0.8 I only ::: 500g
Col_cluster<-c( "#8c96c6" , "#00B3FFFF", "#4d004b" ,  "#238443" , "#f768a1"  , 
                "#88419d"  , "orange"  , "#fde0dd"  , "yellow2" , "#d4b9da",
                "#78c679"  , "#d9f0a3" ,  "#4eb3d3" , "#081d58" ,  "#00FFB2FF",
                "#0868ac" ,  "#fa9fb5" ) # ,  "#FF0000FF", "#00FF19FF"
names(Col_cluster) <- levels(matrix.umap$Cluster)

# I 500g tempora:  9  16 11 5  1  15 10 12 2  3  6  7  13 8  14 4  0 
### ANNOTATION FOR I 500g   ::: reso 0.8
# annotate('text', label = paste0('0'), x =5.87796611978848, y = 4.88927849587492, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13.8323188377105, y = 3.33963703577245, color = 'black', size=8) + annotate('text', label = paste0('2'), x =4.37210505997082, y = -4.76027538448845, color = 'black', size=8) + annotate('text', label = paste0('3'), x =5.70819445999722, y = 10.7381009141621, color = 'black', size=8) + annotate('text', label = paste0('4'), x =12.7964972240007, y = -2.86295538716177, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-16.1321527784097, y = 1.97123647289346, color = 'black', size=8) + annotate('text', label = paste0('6'), x =1.68546183419877, y = -2.72732755013272, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.12566884892418, y = 0.42130648077516, color = 'black', size=8) + annotate('text', label = paste0('8'), x =4.36677349204667, y = -0.11973692435272, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-14.2896370079564, y = -5.95544902753612, color = 'black', size=8) + annotate('text', label = paste0('10'), x =2.55504142903201, y = -12.3286018026813, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-13.8614232970748, y = -1.13772514143639, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.558855462117866, y = -5.97687983631459, color = 'black', size=8) + annotate('text', label = paste0('13'), x =11.92281084236, y = -6.71437844533307, color = 'black', size=8) + annotate('text', label = paste0('14'), x =9.41082696301275, y = 4.3318754344241, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-7.68612899086991, y = 3.75034526249429, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-12.179984067351, y = -3.01200830413311, color = 'black', size=8)

### ANNOTATION FOR I 500g   ::: reso 1
#annotate('text', label = paste0('0'), x =5.99954783615198, y = 4.68034598909802, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13.9337700057022, y = 3.40031792575116, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.90677381476444, y = 10.9051924170025, color = 'black', size=8) + annotate('text', label = paste0('3'), x =12.1844291363024, y = -3.14208903183113, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.07654219662795, y = -2.102784513464, color = 'black', size=8) + annotate('text', label = paste0('5'), x =4.61350702674114, y = -4.30663914020788, color = 'black', size=8) + annotate('text', label = paste0('6'), x =8.1078445857513, y = 0.572450688846357, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-16.4004864236772, y = 1.76309796214115, color = 'black', size=8) + annotate('text', label = paste0('8'), x =4.45782101928132, y = -0.168489255345233, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-14.3147319394345, y = -6.06100308806291, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.9368176922869, y = -1.76549991048023, color = 'black', size=8) + annotate('text', label = paste0('11'), x =2.55164424517314, y = -12.2962653924292, color = 'black', size=8) + annotate('text', label = paste0('12'), x =1.43475197365424, y = -4.50311885404145, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-11.9836467355599, y = 2.53161066011591, color = 'black', size=8) + annotate('text', label = paste0('14'), x =11.8540050443975, y = -6.56503335609062, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-1.36524411838435, y = -8.17723768517781, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.89484791013131, y = 4.69062401508353, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-7.74539474559748, y = 3.72747886425055, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-2.03883463355515, y = -3.55750773737377, color = 'black', size=8)


# I 400g tempora ::
#Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels= c("9", "17", "13", "8" , "0",  "18", "10", "14", "11", "2",  "15", "4",  "12", "6",  "7","3" , "16", "1",  "5"))

### ANNOTATION FOR I 400g   ::: reso 0.8
# annotate('text', label = paste0('0'), x =-14.1918650367615, y = 2.45005014455631, color = 'black', size=8) + annotate('text', label = paste0('1'), x =5.6970107270438, y = 4.12264338026715, color = 'black', size=8) + annotate('text', label = paste0('2'), x =4.17895120936543, y = -6.12718354401615, color = 'black', size=8) + annotate('text', label = paste0('3'), x =3.78523717383112, y = -1.68636092836816, color = 'black', size=8) + annotate('text', label = paste0('4'), x =4.86463958931283, y = 11.402404673131, color = 'black', size=8) + annotate('text', label = paste0('5'), x =12.1072381173036, y = -3.37220108681221, color = 'black', size=8) + annotate('text', label = paste0('6'), x =8.64652473299921, y = 0.894134299920099, color = 'black', size=8) + annotate('text', label = paste0('7'), x =1.7418199645168, y = 3.7102361197927, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-16.0574908073464, y = 1.38469621997143, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-14.5963346180304, y = -5.65967140386441, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-3.67275884416982, y = -8.19837868679402, color = 'black', size=8) + annotate('text', label = paste0('11'), x =2.45427240898523, y = -5.16217934344367, color = 'black', size=8) + annotate('text', label = paste0('12'), x =15.4128077013437, y = -1.06893661370651, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-14.7691349466043, y = -2.24189194560544, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.0627398503428048, y = -10.0340410489269, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-1.94969935167701, y = -4.28180389651811, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.57478942048938, y = 0.135298198103853, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-12.4602405092264, y = -3.70587923633991, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-9.15418392068431, y = 9.64761115250162, color = 'black', size=8)

### ANNOTATION FOR I::: reso 0.8
#annotate('text', label = paste0('0'), x =0.408182061223627, y = 3.0224932614876, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-15.4148570761659, y = 1.09461633809359, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.45550014992074, y = 0.962192760895717, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.633760864213644, y = -0.0364754294229297, color = 'black', size=8) + annotate('text', label = paste0('12'), x =5.85145803293334, y = -1.03049970152282, color = 'black', size=8) + annotate('text', label = paste0('13'), x =11.1984854513172, y = -9.45272637147579, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-4.92915466389885, y = 2.57609514534916, color = 'black', size=8) + annotate('text', label = paste0('15'), x =8.46768359554022, y = 0.865903221023236, color = 'black', size=8) + annotate('text', label = paste0('16'), x =12.5355419082741, y = 6.02785664889035, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-4.58120600251069, y = -16.3291196188194, color = 'black', size=8) + annotate('text', label = paste0('18'), x =6.19793512563697, y = 10.4960847105223, color = 'black', size=8) + annotate('text', label = paste0('2'), x =6.75472576743352, y = -3.01376188489095, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-1, y = 6, color = 'black', size=8) + annotate('text', label = paste0('4'), x =12.0177967751237, y = -4.93485963804348, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.22543471735459, y = 5.5446954312839, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-15.6127975276892, y = 3.5, color = 'black', size=8) + annotate('text', label = paste0('7'), x =4.28978343410821, y = -6.19703549791889, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-3.02900115942241, y = -6.94122543964973, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-12.2561723473624, y = 7.13598374803715, color = 'black', size=8)

### ANNOTATION FOR R::: reso 0.8
# annotate('text', label = paste0('0'), x =-3.45642881523893, y = -1.16865608353008, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9, y = 1.20041811420326, color = 'black', size=8) + annotate('text', label = paste0('2'), x =6.93022842497452, y = 6, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-2.01306296836638, y = -0.148318949858674, color = 'black', size=8) + annotate('text', label = paste0('4'), x =-7.61853096844857, y = 7.37072348119153, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-3, y = -5, color = 'black', size=8) + annotate('text', label = paste0('6'), x =13, y = -2, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.1214536134117, y = 6.0361858560099, color = 'black', size=8) + annotate('text', label = paste0('8'), x =2.77956317997333, y = 4.5, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-12, y = 0.170372736285276, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-10, y = -7, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-1.07578869040965, y = -10.6355224876714, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0280023706045, y = 4.43063015252513, color = 'black', size=8) + annotate('text', label = paste0('13'), x =8.2902770868567, y = -4.12080345312305, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-1.00177421838587, y = -5.88854235521251, color = 'black', size=8) + annotate('text', label = paste0('15'), x =1.70760474758871, y = -5.42889992216071, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.29502755116225, y = -1.14297267439169, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-7.37210218061691, y = -4.23909794401224, color = 'black', size=8) + annotate('text', label = paste0('18'), x =7, y = -7, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-5.5509586359019, y = -12.5376869794186, color = 'black', size=8) + annotate('text', label = paste0('20'), x =11.9769697390055, y = 2.16874421082985, color = 'black', size=8)

# R reso 0.8 500g:::
#annotate('text', label = paste0('0'), x =-1.82181175956801, y = 4.49999935653872, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.85987776275363, y = -4.08730875205322, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-8.77473447623156, y = -2.20602558492533, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-11.9784764388457, y = -3.30252542357977, color = 'black', size=8) + annotate('text', label = paste0('4'), x =12.1904556127542, y = 5.71974924785885, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.75578650613035, y = -0.405585607468961, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-5.59620085133271, y = -4.90009067204302, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-7.92912409998982, y = -6.47019665269102, color = 'black', size=8) + annotate('text', label = paste0('8'), x =12.6922973162172, y = 0.834883693775545, color = 'black', size=8) + annotate('text', label = paste0('9'), x =7.70544237129227, y = 2.44967045106164, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.93339268545303, y = 7.05177242736881, color = 'black', size=8) + annotate('text', label = paste0('11'), x =0.762764834454813, y = 2.82502410907024, color = 'black', size=8) + annotate('text', label = paste0('12'), x =1.66916328116818, y = 5.58051942859595, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-12.2000553204784, y = 2.32987328612614, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-10.1664476877147, y = 0.0955457352512295, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.36852911837374, y = -2.05469576367534, color = 'black', size=8)


# R reso 0.8 400g:::
#annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)



message("----------------------umap.Cluster-----------------------")


umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.3, size=0.2) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("Cluster", values = Col_cluster) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) + 
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.8, size=6))) + 
  theme_classic() +
  #theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)

umap.Cluster

#umap.Cluster_no_legend <- umap.Cluster
#umap.Cluster_w_legend <- umap.Cluster


png(paste0("UMAP_MergedRealigned51_", Project, dataset_name,  "_ggplot_", "umap.Cluster_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Cluster_w_legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project, dataset_name, "_ggplot_", "umap.Cluster_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Cluster_no_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project, dataset_name,  "_ggplot_", "umap.Cluster_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Cluster_w_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project, dataset_name, "_ggplot_", "umap.Cluster_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Cluster_no_legend
dev.off()




message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.1) +
  #geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
  # linetype='dashed', colour="black", bins=3, alpha=0.5) +
  #ggtitle(paste0( dataset_name ,"QC ::: ", reso, " Cell Cycle" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) + coord_fixed() +
  theme_classic() +  #   theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  #theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)



umap.cc


#umap.cc_no_legend <- umap.cc
#umap.cc_w_legend <- umap.cc


png(paste0("UMAP_MergedRealigned51_", Project, dataset_name,  "_ggplot_", "umap.CellCycle_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.cc_w_legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  dataset_name, "_ggplot_", "umap.CellCycle_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.cc_no_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project, dataset_name,  "_ggplot_", "umap.CellCycle_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.cc_w_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project, dataset_name, "_ggplot_", "umap.CellCycle_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.cc_no_legend
dev.off()




message("----------------------umap.Experiment-----------------------")

matrix.umap$Age <- gsub( "*.\\.", "", matrix.umap$Experiment)
unique(matrix.umap$Age)
matrix.umap$Age <- factor(matrix.umap$Age, levels = c("0", "1", "4", "24", "36", "48"))

umap.Experiment    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Age)) +
  geom_point(alpha=0.2, size=0.3) +
  #ggtitle(paste0("Resolution ", reso, " :: Experiment " )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Age, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0"="black", "1"="yellow",  "4"="red", "24"="lightblue2", "36"="blue", "48"="purple"),   limits=c( '0', '1', '4', '24', '36', '48')) + 
  theme_classic() + coord_fixed() +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=6), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) + 
  #theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)



umap.Experiment


#umap.Expt_no_legend <- umap.Experiment
#umap.Expt_w_legend <- umap.Experiment


png(paste0("UMAP_MergedRealigned51_", Project, dataset_name,  "_ggplot_", "umap.Experiment_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Expt_w_legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  dataset_name, "_ggplot_", "umap.Experiment_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Expt_no_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project, dataset_name,  "_ggplot_", "umap.Experiment_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Expt_w_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project, dataset_name, "_ggplot_", "umap.Experiment_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Expt_no_legend
dev.off()








message("----------------------umap.t.0-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Age == "0", "0", "other_timepoints")

table(matrix.umap$Mixer)

umap.t0    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.3, size=0.3) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "0" = "black")) +
  coord_fixed() +
  xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)



 
message("----------------------umap.t.1-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Age == "1", "1", "other_timepoints")

umap.Experiment1    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.2, size=0.2) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "1" = "yellow2")) +
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)





message("----------------------umap.t.4-----------------------")
  
matrix.umap$Mixer <- ifelse(matrix.umap$Age == "4", "4", "other_timepoints")

umap.Experiment4    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.2, size=0.2) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "4" = "red")) +
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)




message("----------------------umap.t.24-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Age == "24", "24", "other_timepoints")

umap.Experiment24    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.2, size=0.2) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "24" = "cyan2")) +
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))  + 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)




   

message("----------------------umap.t.36-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Age == "36", "36", "other_timepoints")

umap.Experiment36    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.2, size=0.2) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "36" = "blue")) +
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)








message("----------------------umap.t.48-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Age == "48", "48", "other_timepoints")

umap.Experiment48    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.2, size=0.2) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "48" = "darkorchid")) +
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =2.6101003654104, y = 5.40826500234414, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.56098292188039, y = -3.68426714162814, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.293478002195, y = -0.706010187561131, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-9.34650142157551, y = -4.69780231730104, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.09545296078564, y = -1.63750734080036, color = 'black', size=8) + annotate('text', label = paste0('5'), x =9.34829102689791, y = -5.97240444870357, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-12.269902783745, y = 2.40351299042826, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-0.893381736135114, y = 3.86210330370172, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-7.35380325368887, y = 7.32036060126618, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-13.5194301697484, y = -1.14238475106098, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-6.53656665537782, y = 2.52958923148266, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.856189159424339, y = 5.00205293355308, color = 'black', size=8) + annotate('text', label = paste0('12'), x =11.2196056845607, y = 3.72038591326236, color = 'black', size=8) + annotate('text', label = paste0('13'), x =10.2169244632044, y = 0.555725676880679, color = 'black', size=8) + annotate('text', label = paste0('14'), x =12.6261373685243, y = 1.6326712497485, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-3.78730236428662, y = 4.39321069958959, color = 'black', size=8) + annotate('text', label = paste0('16'), x =7.590392489613, y = 2.45521931922745, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-16.9077900394019, y = 1.33188122190122, color = 'black', size=8)




pdf(paste0("UMAP_MergedRealigned51_", Project, dataset_name,  "_ggplot_", "umap.Expt_grid_eachTimePoint_","res_", reso,".pdf"), width=20, height=20)
par(bg=NA)
plot_grid(umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 2, align = "hv")
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project, dataset_name, "_ggplot_", "umap.Expt_grid_eachTimePoint_","res_", reso,".png"), width=2000, height=1000,  type="cairo")
par(bg=NA)
plot_grid(umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 2, align = "hv")
dev.off()


list_umaps_for_plotting <- list(umap.t0=umap.t0,
                                umap.Experiment1=umap.Experiment1, 
                                umap.Experiment4=umap.Experiment4,
                                umap.Experiment24=umap.Experiment24,
                                umap.Experiment36=umap.Experiment36, 
                                umap.Experiment48=umap.Experiment48)

saveRDS(list_umaps_for_plotting, paste0("list_umaps_for_plotting_noLegend_", dataset_name, "res0.8.Rds"))

list_umaps_for_plotting2 <- list(umap.cc_no_legend=umap.cc_no_legend, 
                                 umap.cc_w_legend=umap.cc_w_legend, 
                                 umap.Cluster_no_legend=umap.Cluster_no_legend, 
                                 umap.Cluster_w_legend=umap.Cluster_w_legend)

saveRDS(list_umaps_for_plotting2, paste0("list_umaps_for_plotting_2_", dataset_name, "res0.8.Rds"))



pdf(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "9x_GRID_PAPER_FIGURE", reso, ".pdf", sep="_"), width=20,height=18)
#png(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "9x_GRID_PAPER_FIGURE", reso, ".png", sep="_"), width=2500,height=2000, type = "cairo")
par(bg=NA)
plot_grid( list_umaps_for_plotting2[["umap.Cluster_no_legend"]], 
           list_umaps_for_plotting2[["umap.cc_no_legend"]],
           NA, 
           list_umaps_for_plotting[["umap.t0"]], 
           list_umaps_for_plotting[["umap.Experiment1"]],
           list_umaps_for_plotting[["umap.Experiment4"]],
           list_umaps_for_plotting[["umap.Experiment24"]],
           list_umaps_for_plotting[["umap.Experiment36"]],
           list_umaps_for_plotting[["umap.Experiment48"]], ncol = 3, align = "vh")
dev.off()





message("+-------------------------------------------------------------------------------")
message("+                 plot previously identified genes                              ")
message("+-------------------------------------------------------------------------------")


# top QC check genes: 
FeaturePlot(matrix.su_sep,  reduction = "umap", features = c("Elf5", "Eomes", "Sox2","Ascl2","Gcm1","Ovol2", "Anxa1" , "Plac1", "Uba6")) 
FeaturePlot(matrix.su_sep,  reduction = "umap", features = c("Pparg","Gata1", "Klf7","Cdx2","Bhlhe40","Myc", "Bnip3", "Lgals3", "Krt18")) # "Hdac6" , "Sp3", "Tceb2

#FeaturePlot(matrix.su_sep, features = c("Usf2", "Nfya", "Yy1", "Maff", "Gcm1","Ovol2", "Pparg","Ghrl1")) # LP
#FeaturePlot(matrix.su_sep, features = c("Vezf1", "Gata1", "Klf7", "Bdp1","Ascl2")) # EPC
#FeaturePlot(matrix.su_sep, features = c("Flt1", "Anxa1", "Uba6", "Bnip3","Tceb2", "Plac1")) # EPC
#FeaturePlot(matrix.su_sep, features = c("Elf5", "Cdx2","Eomes","Tead4","Sox2","Ascl2","Gcm1","Ovol2"))
# LP ("Gcm1","Ovol2"), EPC = Ascl2FeaturePlot(matrix.su_sep, features = c("Bhlhe40","Myc", "Hdac6","Smarca4",  "Foxo3","Sp3","Srebf2","Tead3", "Ppard")) # TFs









message("--------------------------------------------------------------------------------")
message("+            Find Markers for all cluster pairs res 0.8                         ")
message("+-------------------------------------------------------------------------------")

# switch to RNA assay for finding markers!!!

if(dataset_name == "Remove_29_500g"){
markerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_R_500g_markers"
setwd(markerDir)
  } else if(dataset_name == "Inhibit_32_500g"){
  markerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_I_500g_markers"
  setwd(markerDir)
} else {
  print("no dataset name specified")
  }


# set up variables:
l2fc_cutoff <- 0.25
reso <- "res.0.8"
last_cluster <- length(unique(matrix.su_sep@meta.data[, paste0("SCT_snn_", reso)]))-1

head(matrix.su_sep@meta.data,2)
Idents(matrix.su_sep) <-  matrix.su_sep@meta.data[,paste0("SCT_snn_", reso)]

cluster_numbers <- c(0:(  length(unique(matrix.su_sep@meta.data[, paste0("SCT_snn_", reso)]))-1 ))
marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = length(cluster_numbers), nrow = length(cluster_numbers) ))

marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = last_cluster+1, nrow = last_cluster+1 ))
rownames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15", "c16") #,"c16","c17","c18","c19","c20","c21","c22")
colnames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15", "c16") #,"c16","c17","c18","c19","c20","c21","c22")

table(matrix.su_sep@meta.data[, paste0("SCT_snn_", reso)])



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
                    result_21 = result_21)
names(result_list)

#saveRDS(result_list, "Harmony_R_only_500g_res.0.8__result_list.Rds")




message("--------------------------------------------------------------------------------")
message("+    Now caculate number of DE in pairwise clusters for choser RES              ")
message("+-------------------------------------------------------------------------------")
setwd(ResDir)
result_list <- readRDS("HARMONY_SEPARATE_R_500g_markers/Harmony_R_only_500g_res.0.8__result_list.Rds")
#result_list <- readRDS("HARMONY_SEPARATE_I_500g_markers/Harmony_I_only_500g_res.0.8__result_list.Rds")


l2fc_cutoff <- 0.25
#l2fc_cutoff <- 1


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
  #length(marker_tbl_list[[z]])
}
marker_tbl_list[[1]]


marker_tbl <- as.data.frame(marker_tbl)
length(marker_tbl_list)
for (z in seq_along(marker_tbl_list)){
  marker_tbl[z,] <- c( rep(NA, (last_cluster+1) -length(marker_tbl_list[[z]]) ), marker_tbl_list[[z]])
}
marker_tbl

length(marker_tbl_list[[1]])
length(marker_tbl_list[[2]])
length(marker_tbl_list[[3]])
length(marker_tbl_list[[4]])

max(marker_tbl, na.rm = T) # 79



#httr::set_config(httr::config(ssl_verifypeer = FALSE))

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
      write.csv(tmp_df, paste0("CTR_dscj1_NewAlign51", dataset_name , reso, "_l2fc", l2fc_cutoff,"_Markers_tbl__clusters_", name_y, ".vs.", name_x, ".csv"))
      print(nrow(y[[x]])) }
  }
}

#save_res_tables(x=1, y=result_01)
#save_res_tables(x=4, y=result_01)

GODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_R_500g_GO"
#MarkerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_I_500g_markers"
setwd(GODir)

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

ht3 = Heatmap(as.matrix(marker_tbl),  col = f3, row_title = "", column_title = "Markers between cluster pairs (res 1.25,  abs(l2fc) > 0.25)", show_row_names = TRUE, heatmap_legend_param = list(title = "Number of markers", legend_height = unit(8, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left",  show_heatmap_legend = TRUE) #  width = unit(10, "cm"),
ht3


pdf(paste("Fig__Pairwise_MARKERS", Project, "ComplexHeatmap", dataset_name, reso, "l2fc0.25", ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

#write.csv(marker_tbl, paste0("CTR_dscj1_MergedRealigned51_Pairwise_cluster_marker_tbl__", dataset_name, "_", reso,"_l2fc", l2fc_cutoff,".csv"), quote = FALSE)








message("--------------------------------------------------------------------------------")
message("+                     GO  using enrichR                               ")
message("+-------------------------------------------------------------------------------")


if(dataset_name == "Remove_29_500g"){
  GODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_R_500g_GO"
  setwd(GODir)
} else if(dataset_name == "Inhibit_32_500g"){
  GODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_I_500g_GO"
  setwd(GODir)
} else {
  print("no dataset name specified")
}


library("enrichR")
enrichR_DB <- as.data.frame(listEnrichrDbs())

#l2fc_cutoff <- 0.6
l2fc_cutoff <- 0.25


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
        write.csv(enrichR_RES, paste0(dataset_name, "_reso_", reso, "__", db, "_l2fc", l2fc_cutoff,"_", comparison_name, ".csv"))
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



#### DECIDE ON PARAMETERS:::
#l2fc_cutoff <- 0.6
l2fc_cutoff <- 0.25
database <- "KEGG_2019_Mouse"
# KEGG_2019_Mouse, GO_Biological_Process_2018, GO_Cellular_Component_2018,  GO_Molecular_Function_2018  BioCarta_2016  WikiPathways_2019_Mouse 
database_list <- c("KEGG_2019_Mouse", "GO_Biological_Process_2018", "GO_Cellular_Component_2018",  "GO_Molecular_Function_2018",  "BioCarta_2016",  "WikiPathways_2019_Mouse")
GO_matrix_list <- list()

for (i in seq_along(database_list)){
  database <- database_list[[i]]
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
  GO_matrix_list[[i]] <- GO_matrix
}
names(GO_matrix_list) <- database_list


#GO_matrix <- data.frame(c("xx", "xxx"))
#colnames(GO_matrix) <- paste0(database,".Term")  

#for (z in seq_along(result_list)){
#  res <- result_list[[z]]
#  res_name <- names(result_list)[[z]]
#  for (w in 1:length(res)){
#    GO_matrix <- save_GO_res(w, res, db = database ) 
#  }
#}

#dim(GO_matrix) # 
#GO_matrix[,1] <- gsub(",",";", GO_matrix[,1])

#Kegg_matrix_l2fc0.6 <- GO_matrix
#GOBP_matrix_l2fc0.6  <- GO_matrix
#GOCC_matrix_l2fc0.6  <- GO_matrix
#GOMF_matrix_l2fc0.6  <- GO_matrix
#BioCarta_matrix_l2fc0.6  <- GO_matrix
#WP_matrix_l2fc0.6  <- GO_matrix

Kegg_matrix_l2fc0.25 <- GO_matrix_list[[1]]
GOBP_matrix_l2fc0.25 <- GO_matrix_list[[2]]
GOCC_matrix_l2fc0.25 <- GO_matrix_list[[3]]
GOMF_matrix_l2fc0.25 <- GO_matrix_list[[4]]
BioCarta_matrix_l2fc0.25 <- GO_matrix_list[[5]]
WP_matrix_l2fc0.25 <- GO_matrix_list[[6]]

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

#saveRDS(GO_res_list, paste(dataset_name, "Harmony_Separate", "l2fc_cutoff",l2fc_cutoff, "_GO_res_list.Rds", sep = "_"))
names(GO_res_list)

for(i in seq_along(GO_res_list)){
  GO_res_list[[i]] <- GO_res_list[[i]][-c(1,2),]
}
head(GO_res_list[[1]],2)

for(i in seq_along(GO_res_list)){
  write.csv(GO_res_list[[i]], paste( "NewAlign51_HARMONY_Separate_", dataset_name, reso, names(GO_res_list)[[i]], ".csv", sep = "_"),  quote = FALSE)
}


GOBP_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Remove_29_500g_res.0.8_GOBP_matrix_l2fc0.25_.csv")
GOCC_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Remove_29_500g_res.0.8_GOCC_matrix_l2fc0.25_.csv")
GOMF_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Remove_29_500g_res.0.8_GOMF_matrix_l2fc0.25_.csv")
BioCarta_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Remove_29_500g_res.0.8_BioCarta_matrix_l2fc0.25_.csv")
WP_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Remove_29_500g_res.0.8_WP_matrix_l2fc0.25_.csv")
Kegg_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Remove_29_500g_res.0.8_Kegg_matrix_l2fc0.25_.csv")


GOBP_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Inhibit_32_500g_res.0.8_GOBP_matrix_l2fc0.25_.csv")
GOCC_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Inhibit_32_500g_res.0.8_GOCC_matrix_l2fc0.25_.csv")
GOMF_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Inhibit_32_500g_res.0.8_GOMF_matrix_l2fc0.25_.csv")
BioCarta_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Inhibit_32_500g_res.0.8_BioCarta_matrix_l2fc0.25_.csv")
WP_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Inhibit_32_500g_res.0.8_WP_matrix_l2fc0.25_.csv")
Kegg_matrix_l2fc0.25 <- read.csv("NewAlign51_HARMONY_Separate__Inhibit_32_500g_res.0.8_Kegg_matrix_l2fc0.25_.csv")


GO_res_list <- list(GOBP_matrix_l2fc0.25=GOBP_matrix_l2fc0.25,
                    GOCC_matrix_l2fc0.25=GOCC_matrix_l2fc0.25,
                    GOMF_matrix_l2fc0.25=GOMF_matrix_l2fc0.25,
                    BioCarta_matrix_l2fc0.25=BioCarta_matrix_l2fc0.25,
                    WP_matrix_l2fc0.25=WP_matrix_l2fc0.25,
                    Kegg_matrix_l2fc0.25=Kegg_matrix_l2fc0.25)

for(i in seq_along(GO_res_list)){ print(dim(GO_res_list[[i]])) }


message("--------------------------------------------------------------------------------")
message("+                 which cluster is T0?????                              ")
message("+-------------------------------------------------------------------------------")

head(matrix.su_sep@meta.data,2)

table(matrix.su_sep@meta.data[,c("SCT_snn_res.0.8", "Experiment")])


T0_clust <- 14 # R
#T0_clust <- 9 # I

GO_matrix_for_plotting <- GO_res_list[[5]]
GO_res_name <- names(GO_res_list)[[5]]



for(i in seq_along(GO_res_list)){
#for(i in c(4,5)){
  colnames(GO_res_list[[i]]) <- gsub( "X", "", colnames(GO_res_list[[i]]))
  GO_matrix_for_plotting <- GO_res_list[[i]]
  #rownames(GO_matrix_for_plotting) <- GO_matrix_for_plotting[,1]
  GO_matrix_for_plotting <- GO_matrix_for_plotting[,-1]
  GO_res_name <- names(GO_res_list)[[i]]
  l2fc_cutoff <- as.numeric(as.character(gsub(".*l2fc", "", GO_res_name)))
  #GO_matrix_for_plotting <- GO_matrix_for_plotting[GO_matrix_for_plotting[,1] != c("xx", "xxx"),]

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
  
  if(T0_clust < 10){
    T0_clust_idx1 <- GO_matrix2_clust1 == paste0("0",T0_clust)
  } else { T0_clust_idx1 <- GO_matrix2_clust1 == T0_clust}
  T0_clust_idx2 <- GO_matrix2_clust2 == T0_clust
  T0_columns <- c(colnames(GO_matrix2)[T0_clust_idx2],colnames(GO_matrix2)[T0_clust_idx1])
  
  dim(GO_matrix2)
  
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
  rownames(GO_matrix2) <- gsub( "Homo sapiens.*" ,  "" , rownames(GO_matrix2))

  
  f1 = colorRamp2( c(0, 0.0001, 0.05, 0.5, 1), c("#006d2c",  "#2ca25f", "#66c2a4", "white", "lightgrey"), space = "RGB") 
  ht1 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = T, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = F, cluster_rows = F , row_names_side ="left", width = unit(ncol(GO_matrix2), "cm"),) # width = unit(140, "cm"),
  #ht1
  
  #pdf(paste( Project, "ComplexHeatmap", dataset_name, "Fig__ALL_Pairwise_", GO_res_name, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix2)/2+5), height=nrow(GO_matrix2))
  #par(bg=NA)
  #draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  #dev.off()
  
  dim(GO_matrix2)
  GO_matrix3 <- GO_matrix2[colnames(GO_matrix2) %in% T0_columns]
  dim(GO_matrix3)
  GO_matrix3 <- GO_matrix3[apply(GO_matrix3, 1, function(y) !all(is.na(y))),]
  dim(GO_matrix3)
  GO_matrix3[is.na(GO_matrix3)] <- 1
  
  changeColNames <- as.data.frame(colnames(GO_matrix3))
  colnames(changeColNames) <- "colname1"
  changeColNames$clust1 <- gsub( ".vs.*", "", changeColNames$colname1)
  changeColNames$clust2 <- gsub( ".*vs.", "", changeColNames$colname1)
  changeColNames$clust1 <- as.numeric(as.character(changeColNames$clust1 ))
  changeColNames$clust2 <- as.numeric(as.character(changeColNames$clust2 ))
  changeColNames$new_colname <- ifelse(changeColNames$clust1 < T0_clust, paste0(changeColNames$clust2, ".vs.", changeColNames$clust1) ,  paste0(changeColNames$clust1, ".vs.", changeColNames$clust2))
  
  colnames(GO_matrix3) <- changeColNames$new_colname
  
  ht1 = Heatmap(as.matrix(GO_matrix3),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE ,  row_dend_side = "right", row_names_side ="left", width = unit(ncol(GO_matrix3), "cm")) # width = unit(140, "cm"),
  ht1
  
  pdf(paste( Project, "ComplexHeatmap", dataset_name, "Fig___T0_to_others_", GO_res_name, "l2fc", l2fc_cutoff, "_green.pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix3)/2+12), height=nrow(GO_matrix3)/3)
  par(bg=NA)
  draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
  
}









message("--------------------------------------------------------------------------------")
message("+          Now caculate FIND TOP MARKERS for all comparisons                    ")
message("+-------------------------------------------------------------------------------")

library(biomaRt)
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene_id'), mart = ensembl)          

result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_R_500g_markers/Harmony_R_only_500g_res.0.8__result_list.Rds")
result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_I_500g_markers/Harmony_I_only_500g_res.0.8__result_list.Rds")

l2fc_cutoff <- 0.6
marker_genes <- "genes"

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
      #print(tmp_df$Row.names) 
      return(tmp_df$Row.names)
      }
  }
}

marker_gene_list <- list()
for (z in seq_along(result_list)){
  gene_list <- list()
  res <- result_list[[z]]
  res_name <- names(result_list)[[z]]
  for (w in 1:length(res)){
    gene_list[[w]] <- save_res_tables(w, res)
  }
  marker_gene_list[[z]] <- gene_list
}


#marker_genes_l2fc1.5 <- unique(unlist(unlist(marker_gene_list)))
#marker_genes_l2fc1 <- unique(unlist(unlist(marker_gene_list)))
#marker_genes_l2fc0.6 <- unique(unlist(unlist(marker_gene_list)))
#marker_genes_l2fc0.25 <- unique(unlist(unlist(marker_gene_list)))



mat_clusterMeans <- AverageExpression(matrix.su_sep, assays = "SCT", slot = "scale.data")
mat_clusterMeans_SCT <- mat_clusterMeans[[1]][rownames(mat_clusterMeans[[1]]) %in% marker_genes_l2fc0.6,]
mat_clusterMeans_SCT[1:5, 1:5]
mat_clusterMeans_SCT <- mat_clusterMeans_SCT - rowMeans(mat_clusterMeans_SCT)
min_val <- min(mat_clusterMeans_SCT)
max_val <- max(mat_clusterMeans_SCT)

clust_order <-c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8") # for R 500g res0.8
#clust_order <-c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0") # for I 500g res0.8
mat_clusterMeans_SCT <- mat_clusterMeans_SCT[, match(clust_order, colnames(mat_clusterMeans_SCT))]

# for slot: data:
#f1 = colorRamp2( c(-10, -3, 0, 3, 10), c("#173715", "green3", "grey95", "#bf7acd", "#4b0e81"), space = "LAB") 
# for scale.data:
f1 = colorRamp2( c(min_val, 0, 3, 10), c("green3", "grey95", "#bf7acd", "#4b0e81"), space = "LAB") 
lgd1 = Legend(col_fun = f1, title = "Mean cluster expression", at = c(min_val, max_val )  )

ht1 = Heatmap((mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "Marker genes (abs l2fc > 0.6)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE ,  row_names_side ="left",  row_title_rot = 90, column_title_rot = 0 , na_col = "lightgrey", row_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "bottom", column_names_side = "top", row_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1
ht1

pdf(paste(Project, dataset_name, "Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_clusterOrderedTempora.pdf", sep="_"), onefile=FALSE, width=7, height=45) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


ht1_t = Heatmap(t(mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "Marker genes (abs l2fc > 0.6)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = FALSE ,  row_names_side ="left",  row_title_rot = 90, column_title_rot = 0 ,  na_col = "lightgrey", column_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "top", column_names_side = "bottom",  row_dend_reorder = TRUE, column_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1  
ht1_t

pdf(paste(Project, dataset_name,"t_Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_clusterOrderedTempora2.pdf", sep="_"), onefile=FALSE, width=50, height=7) 
par(bg=NA)
draw(ht1_t, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()








message("--------------------------------------------------------------------------------")
message("+                              feature plots                                    ")
message("+-------------------------------------------------------------------------------")


FeaturePlot(matrix.su_sep,  reduction = "umap", features = c("Elf5", "Eomes", "Sox2","Ascl2","Gcm1","Ovol2", "Anxa1" , "Plac1", "Uba6")) 
FeaturePlot(matrix.su_sep,  reduction = "umap", features = c("Krt18", "Gata3","Jak1","Tead4")) 

FeaturePlot(matrix.su_sep,  reduction = "umap", features = c("Phlda2" ,"Rhox6" , "Rhox9",  "Malat1" ,"Cdkn1c", "Krt18",  "Lgals3", "Rpl22")) 


# diff between cluster 9 and 10
rownames(subset(result_list[[10]][[2]][[1]], abs(avg_logFC) > 0.3))
# diff between cluster 9 and 12
rownames(subset(result_list[[10]][[4]][[1]], abs(avg_logFC) > 0.3))
# diff between cluster 10 and 5
rownames(subset(result_list[[6]][[11]][[1]], abs(avg_logFC) > 0.4))
# diff between cluster 4 and 12
rownames(subset(result_list[[5]][[9]][[1]], abs(avg_logFC) > 0.4))
# diff between cluster 4 and 7
rownames(subset(result_list[[5]][[4]][[1]], abs(avg_logFC) > 0.4))
# diff between cluster 11 and 5
rownames(subset(result_list[[6]][[7]][[1]], abs(avg_logFC) > 0.3))
# diff between cluster 4 and 7
rownames(subset(result_list[[5]][[4]][[1]], abs(avg_logFC) > 0.4))
# diff between cluster 7 and 8
rownames(subset(result_list[[8]][[2]][[1]], abs(avg_logFC) > 0.3))
rownames(subset(result_list[[8]][[4]][[1]], abs(avg_logFC) > 0.4))

FeaturePlot(matrix.su_sep,  reduction = "umap", features = c( "Degs1","Ranbp1", "Wbp5" ,"Bsg", "Gas5", "Ybx1" ,"Atxn10" ,  "Ugp2"  , "Calr") )
# diff between cluster 12 and 4



top_Clust_12 <- c("Sct", "Slc16a1", "Las1l","Tfrc","Krt19","Sin3b","Hbegf","Plac1","Fgfbp1", "Plbd1", "Ctsb", "Calm2", "Igsf8", "Ugp2", "Degs1")
top_Clust_11 <- c("Crip1","Tmsb4x","Lgals3","1600025M17Rik","Cyb5r3","Krt8","Krt18","Anxa2","S100a6", "Cldn4", "Tfpi")
top_Clust_10 <- c( "Peg3","P4hb","Tmem37","Car4","Gjb2","Wfdc2","Igfbp2","F2rl1","Tmem150a", "Fabp3", "P4hb", "Car2")
top_Clust_5 <- c( "Tfrc","Sin3b","Slc2a3","Aldoa","Gjb5")
top_Clust_9 <- c("Rhox6","Rhox9","Uqcrq", "Ybx1")
top_Clust_4 <- c( "Phlda2","Elf3")
top_Clust_7 <- c( "Rhox6", "Dbi", "Atp5j")
top_Clust_other <- c( "Ndrg1", "Tbrg1", "Rbbp7", "Cald1", "Car2", "Calr")

top_markers <- unique(c(top_Clust_12,top_Clust_11,top_Clust_10,top_Clust_5,top_Clust_9,top_Clust_4,top_Clust_7, top_Clust_other))

FeaturePlot(matrix.su_sep,  reduction = "umap", features = top_markers[1:9])
FeaturePlot(matrix.su_sep,  reduction = "umap", features = top_markers[10:18])
FeaturePlot(matrix.su_sep,  reduction = "umap", features = top_markers[19:27])
FeaturePlot(matrix.su_sep,  reduction = "umap", features = top_markers[28:36])
FeaturePlot(matrix.su_sep,  reduction = "umap", features = top_markers[37:45])
FeaturePlot(matrix.su_sep,  reduction = "umap", features = top_markers[46:53])

#Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.0.8
#tmp_tbl <- FindMarkers(matrix.su_sep, ident.1 = "7", only.pos = TRUE)






message("--------------------------------------------------------------------------------")
message("+                dot plots for selected WP pathways                             ")
message("+-------------------------------------------------------------------------------")

WP_genes_df <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/WP_Gene_list_for_dotplots.csv")

WP_ETC <- na.omit(as.character(unique(WP_genes_df$Electron.Transport.Chain.WP295)))
WP_OxP <- as.character(unique(WP_genes_df$Oxidative.phosphorylation.WP1248))
WP_GaG <- as.character(unique(WP_genes_df$Glycolysis.and.Gluconeogenesis.WP157))
WP_CRP <- as.character(unique(WP_genes_df$Cytoplasmic.Ribosomal.Proteins.WP163))

WP_ETC <- WP_ETC[WP_ETC %in% rownames(GetAssayData(matrix.su_sep, assay = "SCT"))]
WP_OxP <- WP_OxP[WP_OxP %in% rownames(GetAssayData(matrix.su_sep, assay = "SCT"))]
WP_GaG <- WP_GaG[WP_GaG %in% rownames(GetAssayData(matrix.su_sep, assay = "SCT"))]
WP_CRP <- WP_CRP[WP_CRP %in% rownames(GetAssayData(matrix.su_sep, assay = "SCT"))]



pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", dataset_name, "WP_ElectronTransportChain",".pdf", sep="_"), width=length(WP_ETC)/5+2, height=6)
par(bg=NA)
DotPlot(matrix.su_sep, features = WP_ETC , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", dataset_name, "WP_OxidativePhosphorylation",".pdf", sep="_"), width=length(WP_OxP)/5+2, height=6)
par(bg=NA)
DotPlot(matrix.su_sep, features = WP_OxP , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", dataset_name, "WP_GlycolysisAndGlucogenesis",".pdf", sep="_"), width=length(WP_GaG)/5+2, height=6)
par(bg=NA)
DotPlot(matrix.su_sep, features = WP_GaG , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", dataset_name, "WP_CytoplasmicRibosomalProteins",".pdf", sep="_"), width=length(WP_CRP)/5+2, height=6)
par(bg=NA)
DotPlot(matrix.su_sep, features = WP_CRP , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



clust_order <-c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0") 
Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.0.8
Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels = clust_order)

DotPlot(matrix.su_sep, features = c("E2f8", "Pparg", "Gcm1", "Gata2", "Msx2","Hif1a","Mef2a",  "Clip1", "Kif1l","Zfp644","Phactr4", "Brms1l", "Ascl2","Fgfbp1", "Gjb3","Phlda2","Gata1","Elf3", "Hes6","Plac1") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






message("--------------------------------------------------------------------------------")
message("+            FindAllMarkers Heatmap!!!                             ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 0.25

Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels= c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0"))

levels(matrix.su_sep@meta.data$SCT_snn_res.0.8)
Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.0.8


ALL_markers <- FindAllMarkers(matrix.su_sep, logfc.threshold = l2fc_cutoff, verbose = TRUE) # "wilcox" == default
ALL_markers <- ALL_markers[order(abs(ALL_markers$avg_logFC), decreasing = TRUE),]

write.csv(ALL_markers, paste0(Project, "_FindAllMarkers_wilcox_l2fc0.25",".csv"))

top_markers <- ALL_markers %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 5) 

pdf(paste(Project, "Heatmap", "top_markers", "tempora.pdf", sep="_"), onefile=FALSE, width=15, height=12) 
par(bg=NA)
DoHeatmap(matrix.su_sep, features = top_markers$gene, group.colors = Col_cluster)
dev.off()


top_markers <- ALL_markers %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 10) 
levels(Idents(matrix.su_sep))

pdf(paste(Project, "Heatmap", "top_10_markers", "tempora.pdf", sep="_"), onefile=FALSE, width=15, height=20) 
par(bg=NA)
DoHeatmap(matrix.su_sep, features = top_markers$gene, group.colors = Col_cluster) # , raster = FALSE
dev.off()











