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
theme_set(theme_cowplot())


baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T0_500g_res0.4"
setwd(baseDir)
Project <- "CTR_dscj1_realign52_T0_het"

#matrix.su_T0 <- readRDS("CTR_dscj1_realign52_matrix.su_T0_HET.Rds")
#matrix.su_T0 <- readRDS("matrix.su_T0_500g_harmony.Rds")
matrix.su_T0 <- readRDS("matrix.su_T0_500g_harmony_3432.Rds")
sampleTable <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/sampleTable_MERGED_53.csv")

#sampleTable <- read.csv("sampleTable_MERGED_53.csv")
unique(matrix.su_T0@meta.data$orig.ident)
Idents(matrix.su_T0) <- matrix.su_T0@meta.data$SCT_snn_res.0.4

DotPlot(matrix.su_T0, features = c( "Junb","Lgals3", "Krt18", "Jund","Rhox6", "Rhox9", "Jun", "Malat1", "Cdkn1c", "Krt8", "Crip1", "Crip2", "Krt7", "Klf6", "Anxa2","Btg2","H2-D1", "Sct", "Lgals1", "Phlda2","Fth1", "Cd63", "Hand1", "Cldn4", "Apoe", "Gm37108", "Wdr89", "Rpl35", "Rpl14-ps1", "Gm9794", "Fkbb3", "Elf5", "Esrrb", "Cdx2", "Eomes", "Lsp1", "Wdr76", "Zfp961" ) , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#"Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1",     "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt"

FeaturePlot(matrix.su_T0, features = c("Rhox6", "Lgals3","Krt8", "Cdkn1c"), ncol = 4, cols = c("#EFEFEF", "blue"), min.cutoff = 'q1')



nrow(matrix.su_T0@meta.data[matrix.su_T0@meta.data$orig.ident == "Batch0004_9617_N721",])
nrow(matrix.su_T0@meta.data[matrix.su_T0@meta.data$orig.ident == "Batch0004_9617_N729",])
min(matrix.su_T0@meta.data$nFeature_RNA)
matrix.su_T0 <- subset(matrix.su_T0, orig.ident != "Batch0004_9617_N729")




matrix.su_T0@meta.data$Experiment <- "0.0"
tmp_tbl <- as.data.frame(table(matrix.su_T0@meta.data[,c("SCT_snn_res.0.4", "Experiment")]))
tmp_tbl <- reshape2::dcast(tmp_tbl, formula = SCT_snn_res.0.4 ~ Experiment)
write.csv(tmp_tbl, "Harmony_T0_500g_TABLE_no_of_cells_per_Cluster_Experiment.csv")




message("--------------------------------------------------------------------------------")
message("+                                Load objects:::                                ")
message("+-------------------------------------------------------------------------------")

matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")
sampleTable <- read.csv("sampleTable_MERGED_53.csv")

matrix.su_T0 <- subset(matrix.su, Experiment == "0.0")
rm(matrix.su)
head(matrix.su_T0@meta.data,2)
matrix.su_T0@meta.data$Phase_All <- matrix.su_T0@meta.data$Phase
matrix.su_T0@meta.data <- matrix.su_T0@meta.data[,c("orig.ident", "Batch", "percent.mt", "nCount_RNA", "nFeature_RNA","Phase_All")]

nrow(matrix.su_T0@meta.data[matrix.su_T0@meta.data$orig.ident2 == "Batch0004_9617_N721",])
nrow(matrix.su_T0@meta.data[matrix.su_T0@meta.data$orig.ident2 == "Batch0004_9617_N729",])
#matrix.su_T0 <- subset(matrix.su_T0, orig.ident2 != "Batch0004_9617_N721")
#matrix.su_T0 <- subset(matrix.su_T0, orig.ident != "Batch0004_9617_N729")

matrix.su_T0_400g <- subset(matrix.su_T0, nFeature_RNA >= 400)
matrix.su_T0_500g <- subset(matrix.su_T0, nFeature_RNA >= 500)



message("--------------------------------------------------------------------------------")
message("+               cell cycle with SEURAT for T.0                                  ")
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

matrix.su_T0 <- CellCycleScoring(matrix.su_T0, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
matrix.su_T0_400g <- CellCycleScoring(matrix.su_T0_400g, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
matrix.su_T0_500g <- CellCycleScoring(matrix.su_T0_500g, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



message("--------------------------------------------------------------------------------")
message("+            Processing with seurat                                             ")
message("+-------------------------------------------------------------------------------")

matrix.su_T0 <- SCTransform(matrix.su_T0, vars.to.regress = "percent.mt", verbose = TRUE)
matrix.su_T0_400g <- SCTransform(matrix.su_T0_400g, vars.to.regress = "percent.mt", verbose = TRUE)
matrix.su_T0_500g <- SCTransform(matrix.su_T0_500g, vars.to.regress = "percent.mt", verbose = TRUE)

#matrix.su_T0 <- FindNeighbors(matrix.su_T0, dims = 1:30)
#matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0, n.start = 100, n.iter = 10, random.seed = 1)  
#matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.1,n.start = 100, n.iter = 10, random.seed = 1) 
#matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.2,n.start = 100, n.iter = 10, random.seed = 1) 
#matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.3,n.start = 100, n.iter = 10, random.seed = 1) 
#matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.4,n.start = 100, n.iter = 10, random.seed = 1) 

matrix.su_T0 <- RunPCA(matrix.su_T0, assay = "SCT", npcs = 30)
matrix.su_T0_400g <- RunPCA(matrix.su_T0_400g, assay = "SCT", npcs = 30)
matrix.su_T0_500g <- RunPCA(matrix.su_T0_500g, assay = "SCT", npcs = 30)

#matrix.su_T0 <- RunUMAP(matrix.su_T0, reduction = "pca", assay = "SCT", dims = 1:30, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 

#saveRDS(matrix.su_T0, "CTR_dscj1_realign52_matrix.su_T0_HET.Rds")
#saveRDS(matrix.su_T0, "CTR_dscj1_realign52_matrix.su_T0_HET_rmBatch0004_9617_N729.Rds")



message("--------------------------------------------------------------------------------")
message("+         HARMONY TO RESCUE BATCH EFFECT                                        ")
message("+-------------------------------------------------------------------------------")

matrix.su_T0 <- RunHarmony(matrix.su_T0, group.by.vars = "Batch", reduction = "pca", assay.use="SCT")
matrix.su_T0 <- RunUMAP(matrix.su_T0, reduction = "harmony", dims = 1:30, min.dist = 0.01)
matrix.su_T0 <- FindNeighbors(matrix.su_T0, reduction = "harmony", dims = 1:30) 
matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0)
matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.1)
matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.2)
matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.3)
matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.4)
matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 0.6)
matrix.su_T0 <- FindClusters(matrix.su_T0, resolution = 1)


matrix.su_T0_400g <- RunHarmony(matrix.su_T0_400g, group.by.vars = "Batch", reduction = "pca", assay.use="SCT")
matrix.su_T0_400g <- RunUMAP(matrix.su_T0_400g, reduction = "harmony", dims = 1:30, min.dist = 0.01)
matrix.su_T0_400g <- FindNeighbors(matrix.su_T0_400g, reduction = "harmony", dims = 1:30) 
matrix.su_T0_400g <- FindClusters(matrix.su_T0_400g, resolution = 0)
matrix.su_T0_400g <- FindClusters(matrix.su_T0_400g, resolution = 0.1)
matrix.su_T0_400g <- FindClusters(matrix.su_T0_400g, resolution = 0.2)
matrix.su_T0_400g <- FindClusters(matrix.su_T0_400g, resolution = 0.3)
matrix.su_T0_400g <- FindClusters(matrix.su_T0_400g, resolution = 0.4)
matrix.su_T0_400g <- FindClusters(matrix.su_T0_400g, resolution = 0.6)
matrix.su_T0_400g <- FindClusters(matrix.su_T0_400g, resolution = 1)

matrix.su_T0_500g <- RunHarmony(matrix.su_T0_500g, group.by.vars = "Batch", reduction = "pca", assay.use="SCT")
matrix.su_T0_500g <- RunUMAP(matrix.su_T0_500g, reduction = "harmony", dims = 1:30, min.dist = 0.01)
matrix.su_T0_500g <- FindNeighbors(matrix.su_T0_500g, reduction = "harmony", dims = 1:30) 
matrix.su_T0_500g <- FindClusters(matrix.su_T0_500g, resolution = 0)
matrix.su_T0_500g <- FindClusters(matrix.su_T0_500g, resolution = 0.1)
matrix.su_T0_500g <- FindClusters(matrix.su_T0_500g, resolution = 0.2)
matrix.su_T0_500g <- FindClusters(matrix.su_T0_500g, resolution = 0.3)
matrix.su_T0_500g <- FindClusters(matrix.su_T0_500g, resolution = 0.4)
matrix.su_T0_500g <- FindClusters(matrix.su_T0_500g, resolution = 0.6)
matrix.su_T0_500g <- FindClusters(matrix.su_T0_500g, resolution = 1)



message("--------------------------------------------------------------------------------")
message("+            After loading matrix.su_T0- do QC                                  ")
message("+-------------------------------------------------------------------------------")

table(matrix.su_T0@meta.data$SCT_snn_res.0.1) # 3
table(matrix.su_T0@meta.data$SCT_snn_res.0.2) # 4
table(matrix.su_T0@meta.data$SCT_snn_res.0.3) # 4
table(matrix.su_T0@meta.data$SCT_snn_res.0.4) # 5

gc(table(matrix.su_T0_400g@meta.data$SCT_snn_res.0.1) # 1
table(matrix.su_T0_400g@meta.data$SCT_snn_res.0.2) # 3
table(matrix.su_T0_400g@meta.data$SCT_snn_res.0.3) # 5
table(matrix.su_T0_400g@meta.data$SCT_snn_res.0.4) # 6 

table(matrix.su_T0_500g@meta.data$SCT_snn_res.0.1) # 1
table(matrix.su_T0_500g@meta.data$SCT_snn_res.0.2) # 3
table(matrix.su_T0_500g@meta.data$SCT_snn_res.0.3) # 4
table(matrix.su_T0_500g@meta.data$SCT_snn_res.0.4) # 5 
#     0    1    2    3    4 
#  2605  613  562  246   28 

table(matrix.su_T0_500g@meta.data$SCT_snn_res.0.4, matrix.su_T0_500g@meta.data$Phase) # 5 



message("--------------------------------------------------------------------------------")
message("+                             Plot UMAP                                         ")
message("+-------------------------------------------------------------------------------")

dataset_name <- "T0.SCT_HARMONY_500g"
seurat_obj <- matrix.su_T0
head(seurat_obj@meta.data,2)

matrix.umap            <- as.data.frame(Embeddings(object=seurat_obj, reduction="umap"))
matrix.umap$Sample_2   <- rownames(matrix.umap) 
matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
rownames(sampleTable)  <- sampleTable$Sample_2 
matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
matrix.umap$Project     <- sampleTable[matrix.umap$Sample_2, ]$Project
head(matrix.umap)

cell_cycle               <- seurat_obj@meta.data %>% dplyr::select((Phase))  # 
colnames(cell_cycle)     <- c("Phase")
matrix.umap$Phase        <- cell_cycle$Phase

reso                <- "res_0.4"
clust               <- seurat_obj@meta.data %>% dplyr::select((SCT_snn_res.0.4))  # 
colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 



message("----------------------umap.Cluster-----------------------")
# c( "#8c96c6" , "#00B3FFFF", "#4d004b" ,  "#238443" , "#f768a1"  , 
#"#88419d"  , "orange"  , "#fde0dd"  , "yellow2" , "#d4b9da",
#"#78c679"  , "#d9f0a3" ,  "#4eb3d3" , "#081d58" ,  "#00FFB2FF",
#"#0868ac" ,  "#fa9fb5" )

umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.4, size=0.3) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("Cluster", values = c("#d4b9da" ,  "magenta2" ,"#00B3FFFF",  "green2" , "mediumblue"  )) +
  ggtitle(paste0(" " )) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() 
umap.Cluster

#umap.Cluster_0.4_300g <- umap.Cluster
#umap.Cluster_0.4_400g <- umap.Cluster
#umap.Cluster_0.4_500g <- umap.Cluster



message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.5) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0( " " )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("Phase", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() 
umap.cc

#umap.cc_0.4_300g <- umap.cc
#umap.cc_0.4_400g <- umap.cc
#umap.cc_0.4_500g <- umap.cc



message("----------------------umap.Batch-----------------------")
   
matrix.umap$Batch <- factor(matrix.umap$Batch, levels = c("A00", "0B0","00C",  "0BC"))

umap.Batch        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Batch)) + geom_point(alpha=0.3, size=0.3) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0(" " )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("Batch", values = c("A00"="black", "0B0"="red", "00C"="green2", "0BC"="blue")) +
  theme(title=element_text(size=10), text=element_text(size=10),  axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.7, size=6))) + theme_classic()
umap.Batch


#umap.Batch_300g <- umap.Batch
#umap.Batch_400g <- umap.Batch
#umap.Batch_500g <- umap.Batch


plot_grid( umap.Cluster_0.4_300g, umap.cc_0.4_300g, umap.Batch_300g, umap.Cluster_0.4_400g, umap.cc_0.4_400g,  umap.Batch_400g, umap.Cluster_0.4_500g, umap.cc_0.4_500g, umap.Batch_500g, labels= c("T0 300g", "", "", "T0 400g", "", "", "T0 500g", "", "" ),  ncol = 3, align = "vh")

#saveRDS(matrix.su_T0, "matrix.su_T0_300g_harmony.Rds")
#saveRDS(matrix.su_T0_400g, "matrix.su_T0_400g_harmony.Rds")
#saveRDS(matrix.su_T0_500g, "matrix.su_T0_500g_harmony.Rds")


pdf(paste0("UMAP_MergedRealigned51_", Project,  "HARMONY_ggplot_", "umap.grid","res_", reso,"3432.pdf"), width=12, height=4)
par(bg=NA)
plot_grid( umap.Cluster,umap.cc, umap.Batch, ncol = 3, align = "vh")
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "HARMONY_ggplot_", "umap.CLuster_CC","res_", reso,"3432.pdf"), width=8, height=4)
par(bg=NA)
plot_grid( umap.Cluster,umap.cc, ncol = 2, align = "vh")
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "HARMONY_ggplot_", "umap.Batch","res_", reso,"3432.pdf"), width=4, height=4)
par(bg=NA)
plot_grid( umap.Batch, ncol = 1, align = "vh")
dev.off()




message("----------------------umap.QC-batch 00C-----------------------")
unique(matrix.umap$Batch) # 0B0 0BC A00 00C

matrix.umap$Sample_2 <- paste(matrix.umap$Batch , "X", matrix.umap$Sample, sep = "_")
matrix.umap$Sample_2 <- gsub( "0B0_X.*", "other", matrix.umap$Sample_2)
matrix.umap$Sample_2 <- gsub( "0BC_X.*", "other", matrix.umap$Sample_2)
matrix.umap$Sample_2 <- gsub( "A00_X.*", "other", matrix.umap$Sample_2)
matrix.umap$Sample_2 <- gsub( "00C_X_", "", matrix.umap$Sample_2)

#"Batch0004_9617_N722" "Batch0004_9617_N726" "Batch0004_9617_N727" "Batch0004_9617_N720" "Batch0004_9617_N719" "Batch0004_9617_N718" "Batch0004_9617_N729"
unique(matrix.umap$Sample_2)
table(matrix.umap[,c("Sample_2", "Cluster" )])

umap.QC    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Sample_2)) +
  geom_point(alpha=0.4, size=0.5) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0(dataset_name, " QC :::  ", reso, " Cluster" )) +
  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("Batch0004_9617_N722"="yellow", "Batch0004_9617_N726"="red", "Batch0004_9617_N727"="green", "Batch0004_9617_N720"="blue", "Batch0004_9617_N719"= "violet", "Batch0004_9617_N729" = "orange", "Batch0004_9617_N718"= "orchid", "other" = "grey")) +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() 
umap.QC


matrix.su_T0@meta.data$Sample <- rownames(matrix.su_T0@meta.data) 
matrix.su_T0@meta.data$Sample <- gsub("_[AGCTN]{12}$","", matrix.su_T0@meta.data$Sample)
matrix.su_T0 <- subset(matrix.su_T0, Sample != "Batch0004_9617_N729")







message("--------------------------------------------------------------------------------")
message("+      Find Markers                     ")
message("+-------------------------------------------------------------------------------")


Idents(matrix.su_T0) <- matrix.su_T0@meta.data$SCT_snn_res.0.4 # 

markers_T0_c0 <- FindMarkers(matrix.su_T0, ident.1 = "0", verbose = TRUE)
markers_T0_c0 <- markers_T0_c0[order(markers_T0_c0$avg_logFC, decreasing = TRUE),]
head(markers_T0_c0,10)

markers_T0_c1 <- FindMarkers(matrix.su_T0, ident.1 = "1", verbose = TRUE)
markers_T0_c1 <- markers_T0_c1[order(markers_T0_c1$avg_logFC, decreasing = TRUE),]
head(markers_T0_c1,10)

markers_T0_c2 <- FindMarkers(matrix.su_T0, ident.1 = "2", verbose = TRUE)
markers_T0_c2 <- markers_T0_c2[order(markers_T0_c2$avg_logFC, decreasing = TRUE),]
head(markers_T0_c2,10)

markers_T0_c3 <- FindMarkers(matrix.su_T0, ident.1 = "3", verbose = TRUE)
markers_T0_c3 <- markers_T0_c2[order(markers_T0_c2$avg_logFC, decreasing = TRUE),]
head(markers_T0_c3,10)

#markers_T0_c0 <- markers_T0_c0[abs(markers_T0_c0$avg_logFC) > 0.6,]
#markers_T0_c1 <- markers_T0_c1[abs(markers_T0_c1$avg_logFC) > 0.6,]
#markers_T0_c2 <- markers_T0_c2[abs(markers_T0_c2$avg_logFC) > 0.6,]

FeaturePlot(matrix.su_T0, features = c("Krt18", "Lgals3", "Malat1", "Rhox6", "Cenpf", "Atrx", "mt_Cytb", "Ubb", "Aldoa"))



message("--------------------------------------------------------------------------------")
message("+            Find Markers for all cluster pairs res 0.4                         ")
message("+-------------------------------------------------------------------------------")

resDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/MARKERS_T0_harmony_res0.4"
#setwd(resDir)


head(matrix.su_T0@meta.data,2)
Idents(matrix.su_T0) <-  matrix.su_T0@meta.data$SCT_snn_res.0.4
cluster_numbers <- c(0:(length(unique(matrix.su_T0@meta.data$SCT_snn_res.0.4))-1))
marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = length(unique(matrix.su_T0@meta.data$SCT_snn_res.0.4)), nrow = length(unique(matrix.su_T0@meta.data$SCT_snn_res.0.4))))

# set up variables:
last_cluster <-length(unique(matrix.su_T0@meta.data$SCT_snn_res.0.4))-1
l2fc_cutoff <- 0.25
reso <- "HARMONY_500g_T0_SCT_snn_res.0.4"

library(purrr)

findMarkers_for_cluster_pair <- function(i){
  if (j == i){
    print("same cluster")
  } else {
    tmp_markers <- FindMarkers(matrix.su_T0, ident.1 = j, ident.2= i, logfc.threshold = l2fc_cutoff, verbose = TRUE)
    print(length(unique(row.names(tmp_markers))))
    return(tmp_markers)
  }
}
safely_findMarkers = safely( findMarkers_for_cluster_pair )


cluster_numbers <- c(0:last_cluster)
j <- "0"
result = map(cluster_numbers, safely_findMarkers)
#saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_00 <- result

cluster_numbers <- c(1:last_cluster)
j <- "1"
result = map(cluster_numbers, safely_findMarkers)
#saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_01 <- result

cluster_numbers <- c(2:last_cluster)
j <- "2"
result = map(cluster_numbers, safely_findMarkers)
#saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_02 <- result

cluster_numbers <- c(3:last_cluster)
j <- "3"
result = map(cluster_numbers, safely_findMarkers)
#saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_03 <- result

cluster_numbers <- c(4:last_cluster)
j <- "4"
result = map(cluster_numbers, safely_findMarkers)
#saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_04 <- result

cluster_numbers <- c(5:last_cluster)
j <- "5"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_05 <- result

cluster_numbers <- c(6:last_cluster)
j <- "6"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_06 <- result


result_list <- list(result_00  = result_00, 
                    result_01 = result_01, 
                    result_02 = result_02, 
                    result_03 = result_03, 
                    result_04 = result_04, 
                    result_05 = result_05,
                    result_06 = result_06)
names(result_list)
#saveRDS(result_list, "result_list_harmony_T0_res.0.4_3432.Rds")

setwd(baseDir)
result_list <- readRDS("result_list_T0_500g_harmony.Rds")



message("--------------------------------------------------------------------------------")
message("+    Now caculate number of DE in pairwise clusters for choser RES              ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 0.6

marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = last_cluster+1, nrow = last_cluster+1 ))
rownames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6")
colnames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6")


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
  marker_tbl[z,] <- c( rep(NA, last_cluster+1-length(marker_tbl_list[[z]]) ), marker_tbl_list[[z]])
}
marker_tbl

max(marker_tbl, na.rm = T) # 10

f2 = colorRamp2( c(0, 1, 2, 3,5), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 

ht2 = Heatmap(as.matrix(marker_tbl),  col = f2, row_title = "", column_title = paste0("Markers between cluster pairs_", reso,  "_absl2fc", l2fc_cutoff), show_row_names = TRUE, heatmap_legend_param = list(title = "Number of markers", legend_height = unit(8, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left",  show_heatmap_legend = TRUE) #  width = unit(10, "cm"),
ht2


pdf(paste("Fig__Pairwise_MARKERS", Project, "ComplexHeatmap",  reso, "_l2fc",l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
par(bg=NA)
draw(ht2, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()



message("--------------------------------------------------------------------------------")
message("+    Now caculate number of DE in pairwise clusters for choser RES              ")
message("+-------------------------------------------------------------------------------")
l2fc_cutoff <- 0.6

get_all_marker_genes  <- function(x,y){
  if (is.data.frame(y[[x]]$result) == TRUE){ 
    tmp_df  <- y[[x]]$result
    tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
    tmp_df <- subset(tmp_df, abs(tmp_df$avg_logFC) > l2fc_cutoff)
    marker_genes <- unique(rownames(tmp_df))
    return(marker_genes)
  }
}


retrieve_markers <- function(x){
  marker_tbl_vals <- list()
  for(i in 1:length(x)){
    marker_genes_from_1_comparison <- get_all_marker_genes(i,x)
    marker_genes_from_1_comparison <- paste(marker_genes_from_1_comparison, sep = "_")
    print(head(marker_genes_from_1_comparison))
    marker_tbl_vals[i] <- as.data.frame(marker_genes_from_1_comparison)
  }
  return(unlist(marker_tbl_vals))
}


marker_genes_list <- list()
for (z in seq_along(result_list)){
  marker_genes_list[[z]] <- retrieve_markers(result_list[[z]])
  length(marker_genes_list[[z]])
}

marker_genes_list <- unlist(marker_genes_list)
length(marker_genes_list)         # 705 // 65
length(unique(marker_genes_list)) # 200 for l2fc0.25  / 24 for l2fc0.6
length(result_list) # 7

#marker_genes_0.6 <- unique(marker_genes_list)
#marker_genes_0.25 <- unique(marker_genes_list)



marker_genes_0.25[marker_genes_0.25 %in% c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt")]
# Krt18  Cdkn1c Hand1  Plet1  Junb  
DotPlot(matrix.su_T0, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
FeaturePlot(matrix.su_T0, features = c("Hand1", "Krt18", "Cdkn1c", "Plet1", "Junb", "Lgals3", "Wdr89", "Rhox9", "Crip1", "Malat1", "Apoe", "Rpl37"))


matrix_for_ht <- GetAssayData(matrix.su_T0, assay = "SCT", slot = "scale.data")
matrix_for_ht[1:5,1:5]
matrix_for_ht <- matrix_for_ht[rownames(matrix_for_ht) %in% marker_genes_0.6,]

split <- as.data.frame(colnames(matrix_for_ht))
colnames(split) <- "cells"
split$cluster <- matrix.su_T0@meta.data[match(split$cells, rownames(matrix.su_T0@meta.data)),]$SCT_snn_res.0.4
matrix_for_ht <- t(as.matrix(matrix_for_ht))
matrix_for_ht <- matrix_for_ht - rowMeans(matrix_for_ht)
max_val <- max(matrix_for_ht[,-ncol(matrix_for_ht)])
min_val <- min(matrix_for_ht[,-ncol(matrix_for_ht)])

f1 = colorRamp2(c(min_val, -1.5,  0, 3, 6,  max_val ) , c("#6000BF","blue", "#000000", "#DB061F", "#F1CD22", "white" ), space = "RGB")
#f3 = colorRamp2( c(0, 1, 2, 5, 10), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 

ht3 = Heatmap(matrix_for_ht,  col = f1, row_title = "", column_title = "T0 500g: Markers between cluster pairs (res 0.4,  abs(l2fc) > 0.6)", show_row_names = FALSE, heatmap_legend_param = list(title = "MeanCentred Expr", legend_height = unit(5, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , split = split$cluster, row_names_side ="left",  show_heatmap_legend = TRUE)
ht3



#pdf(paste("Fig__Pairwise_MARKERS", Project, "ComplexHeatmap",  reso, "_l2fc",l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
pdf(paste("Fig__HARMONY_markers", Project, "ComplexHeatmap",  reso, "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


marker_genes_0.6 <- as.character(marker_genes_0.6)
DefaultAssay(matrix.su_T0) <-"RNA"
FeaturePlot(matrix.su_T0, features = marker_genes_0.6, slot = "data", ncol = 6)

pdf(paste("Fig__HARMONY_markers", "FeaturePlots",  reso, "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=16, height=9)
par(bg=NA)
FeaturePlot(matrix.su_T0, features = marker_genes_0.6, slot = "data", ncol = 6)
dev.off()






timecourse_markers_list <- readRDS( "timecourse_markers_list_harmony.orig500g_l2fc0.25.Rds")

timecourse_markers_T0_R24 <- timecourse_markers_list[[1]]
timecourse_markers_T0_I24 <- timecourse_markers_list[[2]]
timecourse_markers_T0_R48 <- timecourse_markers_list[[3]]
timecourse_markers_T0_I48 <- timecourse_markers_list[[4]]

# rpl, rps - early
# krt late


timecourse_markers_l2fc0.25 <- unique(c(rownames(timecourse_markers_list[[1]]), rownames(timecourse_markers_list[[2]]), rownames(timecourse_markers_list[[3]]), rownames(timecourse_markers_list[[4]])))

timecourse_markers_l2fc0.6 <- unique(c( rownames(timecourse_markers_list[[1]][timecourse_markers_list[[1]]$avg_logFC < -0.6 ,]), 
                                        rownames(timecourse_markers_list[[2]][timecourse_markers_list[[2]]$avg_logFC < -0.6 ,]), 
                                        rownames(timecourse_markers_list[[3]][timecourse_markers_list[[3]]$avg_logFC < -0.6 ,]), 
                                        rownames(timecourse_markers_list[[4]][timecourse_markers_list[[4]]$avg_logFC < -0.6 ,]) ) )

timecourse_markers_l2fc0.25 <- unique(c( rownames(timecourse_markers_list[[1]][timecourse_markers_list[[1]]$avg_logFC < -0.25 ,]), 
                                        rownames(timecourse_markers_list[[2]][timecourse_markers_list[[2]]$avg_logFC < -0.25 ,]), 
                                        rownames(timecourse_markers_list[[3]][timecourse_markers_list[[3]]$avg_logFC < -0.25 ,]), 
                                        rownames(timecourse_markers_list[[4]][timecourse_markers_list[[4]]$avg_logFC < -0.25 ,]) ) )

#marker_genes_0.6 <- unique(marker_genes_list)
#marker_genes_0.25 <- unique(marker_genes_list)

marker_genes_0.25[marker_genes_0.25 %in% timecourse_markers_l2fc0.25]
marker_genes_0.25[marker_genes_0.25 %in% timecourse_markers_l2fc0.6]
marker_genes_0.6[marker_genes_0.6 %in% timecourse_markers_l2fc0.6]

dotplot1 <- DotPlot(matrix.su_T0, features = marker_genes_0.6[marker_genes_0.6 %in% timecourse_markers_l2fc0.6] , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dotplot2 <- DotPlot(matrix.su_T0, features = marker_genes_0.25[marker_genes_0.25 %in% timecourse_markers_l2fc0.6] , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dotplot3 <- DotPlot(matrix.su_T0, features = marker_genes_0.25[marker_genes_0.25 %in% timecourse_markers_l2fc0.25] , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


FeaturePlot(matrix.su_T0, features = c("Hand1", "Krt18", "Cdkn1c", "Plet1", "Junb", "Lgals3", "Wdr89", "Rhox9", "Crip1", "Malat1", "Apoe", "Rpl37"))

umapCoord <- as.data.frame(Embeddings(object = dropseq.integrated[["umap"]]))
umapCoord$UMAP_1 <- umapCoord$UMAP_1 * -1

FeaturePlot(dropseq.integrated, features = c("Lgals3", "Crip1", "Krt18", "Krt8", "Malat1", "Cdkn1c"), ncol = 3)
FeaturePlot(matrix.su_T0, features = c("Lgals3", "Crip1", "Krt18", "Krt8", "Malat1", "Cdkn1c"), ncol = 3)



pdf(paste("Dotplot_Harmony_T0_500g", "res0.4", "T0_markers_l2fc0.6_in_timecourse_markers_l2fc0.6", "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=6, height=6)
par(bg=NA)
dotplot1
dev.off()

pdf(paste("Dotplot_Harmony_T0_500g", "res0.4", "T0_markers_l2fc0.25_in_timecourse_markers_l2fc0.6", "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=12, height=6)
par(bg=NA)
dotplot2
dev.off()

pdf(paste("Dotplot_Harmony_T0_500g", "res0.4", "T0_markers_l2fc0.25_in_timecourse_markers_l2fc0.25", "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=15, height=6)
par(bg=NA)
dotplot3
dev.off()





Idents(matrix.su_T0) <- matrix.su_T0@meta.data$SCT_snn_res.0.4
mat_clusterMeans <- AverageExpression(matrix.su_T0, assays = "SCT", slot = "scale.data")
mat_clusterMeans_SCT <- mat_clusterMeans[[1]][rownames(mat_clusterMeans[[1]]) %in% marker_genes_0.6,]
mat_clusterMeans_SCT[1:4, 1:4]
mat_clusterMeans_SCT <- mat_clusterMeans_SCT - rowMeans(mat_clusterMeans_SCT)
min_val <- min(mat_clusterMeans_SCT)
max_val <- max(mat_clusterMeans_SCT)

# for scale.data:
f1 = colorRamp2( c(-1.5, 0, 1, 2), c("green3", "grey95", "#bf7acd", "#4b0e81"), space = "LAB") 
lgd1 = Legend(col_fun = f1, title = "Mean cluster expression", at = c(min_val, max_val )  )

ht1 = Heatmap((mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "T0 Marker genes (abs l2fc > 0.25)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE ,  row_names_side ="left",  row_title_rot = 0, column_title_rot = 0 , na_col = "lightgrey", row_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "bottom", column_names_side = "top", row_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1
ht1

pdf(paste(Project, "Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_ColClust_3432.pdf", sep="_"), onefile=FALSE, width=4, height=6) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()




ht1_t = Heatmap(t(mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "T0 Marker genes (abs l2fc > 0.25)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE ,  row_names_side ="left",  row_title_rot = 90, column_title_rot = 0 ,  na_col = "lightgrey", column_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "top", column_names_side = "bottom",  row_dend_reorder = TRUE, column_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1  
ht1_t

pdf(paste(Project,"t_Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_ColClust_3432.pdf", sep="_"), onefile=FALSE, width=8, height=6) 
par(bg=NA)
draw(ht1_t, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()




