#!/storage/Software/packages/anaconda3/bin/Rscript

# 
rm(list=ls())
gc()
options(bitmapType='cairo')
gc()


library(cowplot)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g"

Project <- "HARMONY_T24_500g"



#matrix.su_T24 <- readRDS("CTR_dscj1_realign51_matrix.su_T.24_SCT_clust.Rds")
matrix.su_T24 <- readRDS("HARMONY_T24_500g/CTR_dscj1_matrix.su_500g_T.24_SCT_HARMONY_clust.Rds")
sampleTable <- read.csv("sampleTable_MERGED_53.csv")
Idents(matrix.su_T24) <- matrix.su_T24@meta.data$SCT_snn_res.0.4

DotPlot(matrix.su_T24, features= c("Ppargc1a", "Mitf"), assay = "SCT")

tmp_tbl <- as.data.frame(table(matrix.su_T24@meta.data[,c("SCT_snn_res.0.4", "Experiment")]))
tmp_tbl <- reshape2::dcast(tmp_tbl, formula = SCT_snn_res.0.4 ~ Experiment)
write.csv(tmp_tbl, "Harmony_T24_500g_TABLE_no_of_cells_per_Cluster_Experiment.csv")



message("--------------------------------------------------------------------------------")
message("+                                Load objects:::                                ")
message("+-------------------------------------------------------------------------------")

matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")
sampleTable <- read.csv("sampleTable_MERGED_53.csv")
unique(matrix.su@meta.data$Experiment)

nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N721",])
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N729",])
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N729")

matrix.su <- subset(matrix.su, nFeature_RNA >= 500)

matrix.su <- subset(matrix.su, Experiment %in% c("R.24", "I.24"))
head(matrix.su@meta.data,2)
matrix.su@meta.data$Phase_All <- matrix.su@meta.data$Phase
matrix.su@meta.data <- matrix.su@meta.data[,c("orig.ident", "Batch", "percent.mt", "nCount_RNA", "nFeature_RNA","Phase_All", "Experiment", "Treatment")]

setwd(ResDir)

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

matrix.su <- CellCycleScoring(matrix.su, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


message("--------------------------------------------------------------------------------")
message("+            Processing with seurat                                             ")
message("+-------------------------------------------------------------------------------")

unique(matrix.su@meta.data$orig.ident) # 7
table(matrix.su@meta.data[,c("orig.ident", "Treatment")]) # 3 samples R + 4 samples I.


matrix.su <- SCTransform(matrix.su, vars.to.regress = "percent.mt", verbose = TRUE)
#matrix.su <- FindNeighbors(matrix.su, dims = 1:30)
#matrix.su <- FindClusters(matrix.su, resolution = 0, n.start = 100, n.iter = 10, random.seed = 1)  
#matrix.su <- FindClusters(matrix.su, resolution = 0.2, n.start = 100, n.iter = 10, random.seed = 1) 
#matrix.su <- FindClusters(matrix.su, resolution = 0.4, n.start = 100, n.iter = 10, random.seed = 1) 
#matrix.su <- FindClusters(matrix.su, resolution = 0.6, n.start = 100, n.iter = 10, random.seed = 1) 
#matrix.su <- FindClusters(matrix.su, resolution = 0.8, n.start = 100, n.iter = 10, random.seed = 1) 

matrix.su <- RunPCA(matrix.su, assay = "SCT", npcs = 30)
#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT", dims = 1:30, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 

#saveRDS(matrix.su, "CTR_dscj1_realign51_matrix.su_T.24_SCT_clust.Rds")
   



matrix.su_T24 <- FindClusters(matrix.su_T24, resolution = 0.1, n.start = 100, n.iter = 10, random.seed = 1) 

message("--------------------------------------------------------------------------------")
message("+                        POST procesing QC                                      ")
message("+-------------------------------------------------------------------------------")



message("--------------------------------------------------------------------------------")
message("+                  DotPlot Myriams markers                                      ")
message("+-------------------------------------------------------------------------------")

reso <- "SCT_snn_res.0.4"
Idents(matrix.su_T24) <- matrix.su_T24@meta.data[reso]

# ALL sig
DotPlot(matrix.su_T24, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",     "Prl3d1","Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2",  "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Junb", "Fosl1", "Arnt") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# "Prl8a9", "Prl2c2", "Tpbpa","Syna", 

pdf(paste("DotPlot_", Project, "_TSC_diff_markers_", "res_", reso,".pdf", sep = "_"), width=9, height=5)
par(bg=NA)
DotPlot(matrix.su_T24, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",     "Prl3d1","Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2",  "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Junb", "Fosl1", "Arnt") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


FeaturePlot(matrix.su_T24, features = c("Hand1", "Gjb3", "Ascl2", "Plac1", "Ly6e", "Gcm1", "Junb", "Cdkn1c"), ncol = 4)


head(matrix.su_T24@meta.data,2)
tmp_df <- as.data.frame((matrix.su_T24@meta.data[,c("Treatment", "SCT_snn_res.0.4")]))
tmp_df <- subset(tmp_df, tmp_df$Treatment != "0")

tmp_df$Treatment <- as.factor(tmp_df$Treatment)
ggplot(data=tmp_df, aes(SCT_snn_res.0.4)) + geom_bar(aes(fill= (Treatment),y = (..count..)/sum(..count..)), position="dodge")




message("--------------------------------------------------------------------------------")
message("+                        HARMONY - optional                                     ")
message("+-------------------------------------------------------------------------------")
matrix.su_T24 <- matrix.su
table(matrix.su_T24@meta.data$Batch)
library(harmony)
matrix.su_T24 <- RunHarmony(matrix.su_T24, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", reference_values = "ABC")
matrix.su_T24 <- RunUMAP(matrix.su_T24, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors= 30L, spread = 2L, repulsion.strength = 2)

matrix.su_T24 <- FindNeighbors(matrix.su_T24, reduction = "harmony", dims = 1:30) 
matrix.su_T24 <- FindClusters(matrix.su_T24, resolution = 0.2)
matrix.su_T24 <- FindClusters(matrix.su_T24, resolution = 0.4)
matrix.su_T24 <- FindClusters(matrix.su_T24, resolution = 0.6)
matrix.su_T24 <- FindClusters(matrix.su_T24, resolution = 1)
matrix.su_T24 <- FindClusters(matrix.su_T24, resolution = 0.1)

#saveRDS(matrix.su_T24, "CTR_dscj1_matrix.su_500g_T.24_SCT_HARMONY_clust.Rds")


message("--------------------------------------------------------------------------------")
message("+            After loading matrix.su_T24- do QC                                  ")
message("+-------------------------------------------------------------------------------")

unique(matrix.su_T24@meta.data$SCT_snn_res.0)   # 1
unique(matrix.su_T24@meta.data$SCT_snn_res.0.2) # 8
unique(matrix.su_T24@meta.data$SCT_snn_res.0.4) # 12
unique(matrix.su_T24@meta.data$SCT_snn_res.0.1) # 
head(matrix.su_T24@meta.data,2)


message("+-------------------------------------------------------------------------------")
message("+                                 Clustree                                      ")
message("+-------------------------------------------------------------------------------")

tree_to_plot <- clustree::clustree(matrix.su_T24, prefix = "SCT_snn_res.") # "RNA_snn_res."

pdf(paste("matrix.su_REALIGNED_Harmony_T24_500g_", "___clustree.pdf", sep=""), width=12,height=8)
par(bg=NA)
tree_to_plot
dev.off()


message("--------------------------------------------------------------------------------")
message("+                             Plot UMAP                                         ")
message("+-------------------------------------------------------------------------------")

dataset_name <- "T24.SCT_HARMONY_500g"
#dataset_name <- "T24.SCT"
head(matrix.su_T24@meta.data,2)

matrix.umap            <- as.data.frame(Embeddings(object=matrix.su_T24, reduction="umap"))
matrix.umap$Sample_2   <- rownames(matrix.umap) 
matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
rownames(sampleTable)  <- sampleTable$Sample_2 
matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
matrix.umap$Project     <- sampleTable[matrix.umap$Sample_2, ]$Project
head(matrix.umap)

cell_cycle               <- matrix.su_T24@meta.data %>% dplyr::select((Phase))  # 
colnames(cell_cycle)     <- c("Phase")
matrix.umap$Phase        <- cell_cycle$Phase

reso                <- "res_0.4"
clust               <- matrix.su_T24@meta.data %>% dplyr::select((SCT_snn_res.0.4))  # 
colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 



message("----------------------umap.Cluster-----------------------")

umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.4, size=0.3) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0( " " )) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=FALSE, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() 
umap.Cluster

#umap.Cluster_0.2 <- umap.Cluster
#umap.Cluster_0.4 <- umap.Cluster
#umap.Cluster_0.6 <- umap.Cluster
#umap.Cluster_0.1 <- umap.Cluster

plot_grid(umap.Cluster_0.2,umap.Cluster_0.4, umap.Cluster_0.6, align = "hv", ncol = 3)




message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.5) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0( " " )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("Phase", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +  guides(colour = guide_legend(reverse=TRUE, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() 
umap.cc



message("----------------------umap.Batch-----------------------")
unique(matrix.umap$Batch)
matrix.umap$Batch <- factor(matrix.umap$Batch, levels = c("AB0", "A00",  "ABC"))

umap.Batch        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Batch)) + geom_point(alpha=0.3, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0(" " )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("Batch", values = c( "ABC"="blue","A00"="green", "AB0"="red")) +
  theme(title=element_text(size=10), text=element_text(size=10),  axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + theme_classic()
umap.Batch



message("----------------------umap.Experiment-----------------------")

umap.Experiment    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.4, size=0.3) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0( " " )) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=TRUE, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() 
umap.Experiment




pdf(paste0("UMAP_MrgdRealign51_", Project,  "HARMONY_ggplot_", "umap.grid","res_", reso,".pdf"), width=12, height=8)
par(bg=NA)
plot_grid(umap.Cluster_0.4, umap.cc, umap.Batch, umap.Experiment, align = "hv")
dev.off()






message("--------------------------------------------------------------------------------")
message("+            Find Markers for all cluster pairs res 0.4                        ")
message("+-------------------------------------------------------------------------------")

MarkerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g/MARKERS_T24_500g_res0.4"
setwd(MarkerDir)
reso <- "SCT_snn_res.0.4"
  
head(matrix.su_T24@meta.data,2)
Idents(matrix.su_T24) <-  matrix.su_T24@meta.data[[reso]]
cluster_numbers <- c(0:(length(unique(matrix.su_T24@meta.data[[reso]]))-1))
marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = length(unique(matrix.su_T24@meta.data[[reso]])), nrow = length(unique(matrix.su_T24@meta.data[[reso]]))))

# set up variables:
last_cluster <-length(unique(matrix.su_T24@meta.data[[reso]]))-1
l2fc_cutoff <- 0.25
reso <- "T24_500g_SCT_snn_res.0.4"

library(purrr)

findMarkers_for_cluster_pair <- function(i){
  if (j == i){
    print("same cluster")
  } else {
    tmp_markers <- FindMarkers(matrix.su_T24, ident.1 = j, ident.2= i, logfc.threshold = l2fc_cutoff, verbose = TRUE)
    print(length(unique(row.names(tmp_markers))))
    return(tmp_markers)
  }
}
safely_findMarkers = safely( findMarkers_for_cluster_pair )


cluster_numbers <- c(0:last_cluster)
j <- "0"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_00 <- result

cluster_numbers <- c(1:last_cluster)
j <- "1"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_01 <- result

cluster_numbers <- c(2:last_cluster)
j <- "2"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_02 <- result

cluster_numbers <- c(3:last_cluster)
j <- "3"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_03 <- result

cluster_numbers <- c(4:last_cluster)
j <- "4"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
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

cluster_numbers <- c(7:last_cluster)
j <- "7"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_07 <- result

cluster_numbers <- c(8:last_cluster)
j <- "8"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_08 <- result

cluster_numbers <- c(9:last_cluster)
j <- "9"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_09 <- result


result_list <- list(result_00 = result_00, 
                    result_01 = result_01, 
                    result_02 = result_02, 
                    result_03 = result_03, 
                    result_04 = result_04, 
                    result_05 = result_05,
                    result_06 = result_06,
                    result_07 = result_07,
                    result_08 = result_08,
                    result_09 = result_09)
names(result_list)
#saveRDS(result_list, "result_list_HARMONY_T24_500g_SCT_res.0.4.Rds")


#Idents(matrix.su_T24) <- matrix.su_T24@meta.data$Treatment
#I_vs_R_markers <- FindMarkers(matrix.su_T24, ident.1 = "I", ident.2= "R", logfc.threshold = l2fc_cutoff, verbose = TRUE)
#Idents(matrix.su_T24) <- matrix.su_T24@meta.data$SCT_snn_res.0.2

FeaturePlot(matrix.su_T24, features = rownames(I_vs_R_markers))

result_list <- readRDS( "result_list_HARMONY_T24_500g_SCT_res.0.4.Rds")



message("--------------------------------------------------------------------------------")
message("+    Now caculate number of DE in pairwise clusters for choser RES              ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 0.6

marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = last_cluster+1, nrow = last_cluster+1 ))
rownames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6", "c7", "c8", "c9")
colnames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6", "c7", "c8", "c9")
cluster_numbers <- c(0:(length(unique(matrix.su_T24@meta.data$SCT_snn_res.0.4))-1))


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
  marker_tbl[z,] <- c( rep(NA, length(cluster_numbers)-length(marker_tbl_list[[z]]) ), marker_tbl_list[[z]])
}
marker_tbl

max(marker_tbl, na.rm = T) # 27

f2 = colorRamp2( c(0, 1, 2, 3,5), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 

ht2 = Heatmap(as.matrix(marker_tbl),  col = f2, row_title = "", column_title = paste0("Markers between cluster pairs (Res 0.2)"," absl2fc", l2fc_cutoff), show_row_names = TRUE, heatmap_legend_param = list(title = "Number of markers", legend_height = unit(8, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left",  show_heatmap_legend = TRUE)
ht2

reso <- 0.4
pdf(paste("Fig_T24_500g_HARMONY_Pairwise_MARKERS_ComplexHeatmap_",  reso, "_l2fc",l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=6, height=5)
par(bg=NA)
draw(ht2, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()





message("--------------------------------------------------------------------------------")
message("+    Now caculate number of DE in pairwise clusters for choser RES              ")
message("+-------------------------------------------------------------------------------")
result_list <- readRDS(paste0(ResDir,"/MARKERS_T24_500g_res0.4/result_list_HARMONY_T24_500g_SCT_res.0.4.Rds"))

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
length(marker_genes_list)         # 2090 // 2923
length(unique(marker_genes_list)) # 407 for l2fc0.25  / 78 for l2fc0.6
length(result_list) # 8

#marker_genes_0.6 <- unique(marker_genes_list)
#marker_genes_0.25 <- unique(marker_genes_list)

pdf(paste("Fig_T24__HARMONY_DotPlot_markers_ComplexHeatmap_SCT_snn_res", reso, "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=20, height=6)
par(bg=NA)
DotPlot(matrix.su_T24, features = marker_genes_0.6, cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




mat_clusterMeans <- AverageExpression(matrix.su_T24, assays = "SCT", slot = "scale.data")
mat_clusterMeans_SCT <- mat_clusterMeans[[1]][rownames(mat_clusterMeans[[1]]) %in% marker_genes_0.6,]
mat_clusterMeans_SCT[1:5, 1:5]
mat_clusterMeans_SCT <- mat_clusterMeans_SCT - rowMeans(mat_clusterMeans_SCT)
min_val <- min(mat_clusterMeans_SCT)
max_val <- max(mat_clusterMeans_SCT)

f1 = colorRamp2(c(min_val,  0, 3, 6) , c("blue", "#000000", "#DB061F", "#F1CD22" ), space = "RGB")
#f3 = colorRamp2( c(0, 1, 2, 5, 10), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 

ht3 = Heatmap(t(mat_clusterMeans_SCT),  col = f1, row_title = "", column_title = "T24 500g: Markers between cluster pairs (res 0.4,  abs(l2fc) > 0.6), cluster means shown", show_row_names = TRUE, heatmap_legend_param = list(title = "MeanCentred Expr", legend_height = unit(5, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , row_names_side ="left",  show_heatmap_legend = TRUE)
ht3


pdf(paste("Fig_T24__HARMONY_ClusterMeans_markers_ComplexHeatmap_SCT_snn_res", reso, "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=15, height=10)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()





matrix_for_ht <- GetAssayData(matrix.su_T24, assay = "SCT", slot = "scale.data")
matrix_for_ht[1:5,1:5]
matrix_for_ht <- matrix_for_ht[rownames(matrix_for_ht) %in% marker_genes_0.6,]
matrix_for_ht <- as.data.frame(matrix_for_ht)

split <- as.data.frame(colnames(matrix_for_ht))
colnames(split) <- "cells"
split$cluster <- matrix.su_T24@meta.data[match(split$cells, rownames(matrix.su_T24@meta.data)),]$SCT_snn_res.0.4
matrix_for_ht <- t(as.matrix(matrix_for_ht))
matrix_for_ht <- matrix_for_ht - rowMeans(matrix_for_ht)

f1 = colorRamp2(c(min(matrix_for_ht[,-ncol(matrix_for_ht)]), -1.5,  0, 3, 6,   max(matrix_for_ht[,-ncol(matrix_for_ht)])) , c("#6000BF","blue", "#000000", "#DB061F", "#F1CD22", "white" ), space = "RGB")
f3 = colorRamp2( c(0, 1, 2, 5, 10), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 

ht3 = Heatmap(matrix_for_ht,  col = f1, row_title = "", column_title = "T24: Markers between cluster pairs (res 0.2,  abs(l2fc) > 0.6)", show_row_names = FALSE, heatmap_legend_param = list(title = "MeanCentred Expr", legend_height = unit(5, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , split = split$cluster, row_names_side ="left",  show_heatmap_legend = TRUE)
#ht3

dim(matrix_for_ht) # 45756    54

pdf(paste("Fig_T24__HARMONY_markers_ComplexHeatmap_SCT_snn_res", reso, "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=15, height=10)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()



marker_genes_0.6 <- as.character(marker_genes_0.6)
length(marker_genes_0.6) # 60 !
#FeaturePlot(matrix.su_T0, features = marker_genes_0.6, slot = "data", ncol = 6)

pdf(paste("Fig_T24__HARMONY_markers_FeaturePlots_", reso, "_l2fc", l2fc_cutoff, "_set1.pdf", sep="_"), onefile=FALSE, width=16, height=9)
par(bg=NA)
FeaturePlot(matrix.su_T24, features = marker_genes_0.6[1:20],  ncol = 5)
dev.off()

pdf(paste("Fig_T24__HARMONY_markers_FeaturePlots_", reso, "_l2fc", l2fc_cutoff, "_set2.pdf", sep="_"), onefile=FALSE, width=16, height=9)
par(bg=NA)
FeaturePlot(matrix.su_T24, features = marker_genes_0.6[21:40], slot = "data", ncol = 5)
dev.off()

pdf(paste("Fig_T24__HARMONY_markers_FeaturePlots_", reso, "_l2fc", l2fc_cutoff, "_set3.pdf", sep="_"), onefile=FALSE, width=16, height=9)
par(bg=NA)
FeaturePlot(matrix.su_T24, features = marker_genes_0.6[41:60], slot = "data", ncol = 5)
dev.off()

pdf(paste("Fig_T24__HARMONY_markers_FeaturePlots_", reso, "_l2fc", l2fc_cutoff, "_set4.pdf", sep="_"), onefile=FALSE, width=16, height=9)
par(bg=NA)
FeaturePlot(matrix.su_T24, features = marker_genes_0.6[61:78], slot = "data", ncol = 5)
dev.off()









message("--------------------------------------------------------------------------------")
message("+    Now caculate number of DE in pairwise clusters for choser RES              ")
message("+-------------------------------------------------------------------------------")

get_res_tables <- function(x,y){
  if (is.data.frame(y[[x]]$result) == FALSE){ print(0)
  } else {
    y[[x]] <- y[[x]]$result
    tmp_df <- y[[x]]
    tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
    number_to_add <- length(unique(rownames(tmp_df)))
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
  marker_tbl[z,] <- c( rep(NA, 8-length(marker_tbl_list[[z]]) ), marker_tbl_list[[z]])
}
marker_tbl

max(marker_tbl, na.rm = T) # 79


library(biomaRt)
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene_id'), mart = ensembl)          




all_marker_genes_list <- list()

save_res_tables <- function(x,y){
  name_y <- res_name
  name_y <- gsub( "result_","" , name_y)
  name_x <- (x - 1 + as.numeric(name_y))
  if (is.data.frame(y[[x]]$result) == FALSE){ print(0)
  } else {
    #y[[x]] <- y[[x]]$result
    tmp_df <- y[[x]]$result
    if(tmp_df == "same cluster"){  print(0)
    } else { 
      tmp_df <- unique(merge(tmp_df, ensEMBL2id, by.x = "row.names", by.y = "external_gene_name", all.x = TRUE))
      tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
      markers <- tmp_df$Row.names
      print(markers)
      write.csv(tmp_df, paste0("CTR_dscj1_NewAlign52_T24_500g_res" , reso, "Markers_tbl__clusters_", name_y, ".vs.", name_x, "l2fc0.25.csv"))
      print(nrow(tmp_df)) 
      #all_marker_genes_list <- append(all_marker_genes_list, markers)
      return(markers)
      }
  }
}

#save_res_tables(x=1, y=result_01)
#save_res_tables(x=4, y=result_01)

all_marker_genes_list2 <- list()
for (z in seq_along(result_list)){
  res <- result_list[[z]]
  res_name <- names(result_list)[[z]]
  for (w in 1:length(res)){
    all_marker_genes_list[[w]] <- save_res_tables(w, res)
  }
  all_marker_genes_list2[[z]] <- all_marker_genes_list
}

all_marker_genes_list2 <- unique(unlist(all_marker_genes_list))
all_marker_genes_list2 <- all_marker_genes_list2[all_marker_genes_list2 != "0"]
#marker_genes_0.6 <- all_marker_genes_list2
#marker_genes_0.25 <- all_marker_genes_list2


library(ComplexHeatmap)
library(circlize)


f3 = colorRamp2( c(0, 1, 2, 5, 10), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 
f3 = colorRamp2( c(min(matrix_for_ht), -2,  mean(matrix_for_ht), 10,max(matrix_for_ht)), c("olivedrab","green","white","darkorchid", "darkorchid4" ), space = "LAB") 

ht3 = Heatmap(as.matrix(marker_tbl),  col = f3, row_title = "", column_title = "Markers between cluster pairs (res 0.4,  abs(l2fc) > 0.6)", show_row_names = TRUE, heatmap_legend_param = list(title = "Number of markers", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left",  show_heatmap_legend = TRUE) #  width = unit(10, "cm"),
ht3


pdf(paste("Fig__Pairwise_MARKERS", Project, "ComplexHeatmap",  reso, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

#write.csv(marker_tbl, "CTR_dscj1_MergedRealigned52_30L_MARKERS___Pairwise_cluster_marker_tbl__res.1_l2f0.6.csv", quote = FALSE)



matrix_for_ht <- GetAssayData(matrix.su_T24, assay = "SCT", slot = "scale.data")
matrix_for_ht[1:5,1:5]
matrix_for_ht <- matrix_for_ht[rownames(matrix_for_ht) %in% marker_genes_0.6,]

split <- as.data.frame(colnames(matrix_for_ht))
colnames(split) <- "cells"
split$cluster <- matrix.su_T24@meta.data[match(split$cells, rownames(matrix.su_T24@meta.data)),]$SCT_snn_res.0.4
matrix_for_ht <- t(as.matrix(matrix_for_ht))
matrix_for_ht <- matrix_for_ht - rowMeans(matrix_for_ht)

f3 = colorRamp2( c(min(matrix_for_ht), -2,  mean(matrix_for_ht), 10,max(matrix_for_ht)), c("olivedrab","green","white","darkorchid", "darkorchid4" ), space = "LAB") 


ht3 = Heatmap(matrix_for_ht,  col = f3, row_title = "", column_title = "T24: Markers between cluster pairs (res 0.4,  abs(l2fc) > 0.6)", show_row_names = FALSE, heatmap_legend_param = list(title = "MeanCentred Expr", legend_height = unit(5, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , split = split$cluster, row_names_side ="left",  show_heatmap_legend = TRUE)
ht3

#ht3_0.6 <- ht3

png(paste( Project, "ComplexHeatmap", "_Pairwise_Markers_",  "res.0.4", "l2fc", l2fc_cutoff, ".png", sep="_"), width=1000, height=2000, type = "cairo")
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


FeaturePlot(matrix.su_T24, features = c("Saa3", "Ccl2", "Igf2", "Ybx1", "Sparc", "Atp5b"), ncol = 4)












message("--------------------------------------------------------------------------------")
message("+                     GO  using enrichR                               ")
message("+-------------------------------------------------------------------------------")

markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g/GO_MARKERS_T24_500g_res0.4"
setwd(markerGODir)

library("enrichR")
enrichR_DB <- as.data.frame(listEnrichrDbs())

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
    #tmp_df <- unique(merge(tmp_df, ensEMBL2id, by.x = "row.names", by.y = "external_gene_name", all.x = TRUE))
    tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
    tmp_df <- subset(tmp_df, abs(tmp_df[,2]) > l2fc_cutoff)
    if( nrow(tmp_df) > 1){  
      enrichR_RES <- as.data.frame(enrichr(rownames(tmp_df), databases = db))
      enrichR_RES <- subset(enrichR_RES, enrichR_RES[,4] < 0.05)
      print(nrow(enrichR_RES))
      if( nrow(enrichR_RES) > 0 ){
        write.csv(enrichR_RES, paste0("NewAlign_T24_500g_SCT__Res__0.2__CCregr_", db, "_l2fc", l2fc_cutoff,"_", comparison_name, ".csv"))
        GO_matrix_to_add <- enrichR_RES[,c(1,4)]
        colnames(GO_matrix_to_add)[2] <- comparison_name
        GO_matrix <- unique(merge(GO_matrix, GO_matrix_to_add, by.x = db.term, by.y =db.term, all.x = TRUE, all.y = TRUE ))
        print(head(GO_matrix,3))
      } 
    }
  }
  return(GO_matrix)
}


#### DECIDE ON PARAMETERS:::
#l2fc_cutoff <- 0.6
l2fc_cutoff <- 0.25
database <- "KEGG_2019_Mouse"
# GO_Biological_Process_2018, GO_Cellular_Component_2018,  GO_Molecular_Function_2018  BioCarta_2016  WikiPathways_2019_Mouse KEGG_2019_Mouse

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

#GOBP_matrix_l2fc0.25 <- GO_matrix
#GOCC_matrix_l2fc0.25 <- GO_matrix
#GOMF_matrix_l2fc0.25 <- GO_matrix
#BioCarta_matrix_l2fc0.25 <- GO_matrix
#WP_matrix_l2fc0.25 <- GO_matrix
#Kegg_matrix_l2fc0.25 <- GO_matrix

#GOBP_matrix_l2fc0.6  <- GO_matrix
#GOCC_matrix_l2fc0.6  <- GO_matrix
#GOMF_matrix_l2fc0.6  <- GO_matrix
#BioCarta_matrix_l2fc0.6  <- GO_matrix
#WP_matrix_l2fc0.6  <- GO_matrix
#Kegg_matrix_l2fc0.6 <- GO_matrix

GO_res_list <- list(GOBP_matrix_l2fc0.25=GOBP_matrix_l2fc0.25,
                    GOCC_matrix_l2fc0.25=GOCC_matrix_l2fc0.25,
                    GOMF_matrix_l2fc0.25=GOMF_matrix_l2fc0.25,
                    BioCarta_matrix_l2fc0.25=BioCarta_matrix_l2fc0.25,
                    WP_matrix_l2fc0.25=WP_matrix_l2fc0.25,
                    Kegg_matrix_l2fc0.25=Kegg_matrix_l2fc0.25)

#saveRDS(GO_res_list, "HARMONY_T24_500g_SCT_res0.4_GO_res_list_l2fc0.25.Rds")
#GO_res_list <- readRDS( "HARMONY_T24_500g_SCT_res0.4_GO_res_list_l2fc0.25.Rds")

names(GO_res_list)

GO_res_list <- list(GOBP_matrix_l2fc0.6=GOBP_matrix_l2fc0.6,
                    GOCC_matrix_l2fc0.6=GOCC_matrix_l2fc0.6,
                    GOMF_matrix_l2fc0.6=GOMF_matrix_l2fc0.6,
                    BioCarta_matrix_l2fc0.6=BioCarta_matrix_l2fc0.6,
                    WP_matrix_l2fc0.6=WP_matrix_l2fc0.6,
                    Kegg_matrix_l2fc_0.6=Kegg_matrix_l2fc_0.6)


for(i in seq_along(GO_res_list)){
  GO_res_list[[i]] <- GO_res_list[[i]][-c(1,2),]
}
head(GO_res_list[[1]],2)

for(i in seq_along(GO_res_list)){
  write.csv(GO_res_list[[i]], paste( "Harmony_T24_500g_SCT_res_0.4_l2fc0.25", names(GO_res_list)[[i]], ".csv", sep = "_"),  quote = FALSE)
}




for(i in seq_along(GO_res_list)){
#for(i in c(4,5)){
  GO_matrix_for_plotting <- GO_res_list[[i]]
  GO_res_name <- names(GO_res_list)[[i]]
  l2fc_cutoff <- as.numeric(as.character(gsub(".*l2fc", "", GO_res_name)))

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


  f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
  ht1 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = T, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = F, cluster_rows = F , row_names_side ="left", width = unit(ncol(GO_matrix2), "cm"),) # width = unit(140, "cm"),
  #ht1

  pdf(paste( Project, "Harmony_ComplexHeatmap", "Fig__ALL_Pairwise_", GO_res_name,  "res.0.4", "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix2)+5), height=nrow(GO_matrix2)/3)
  par(bg=NA)
  draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
}







message("--------------------------------------------------------------------------------")
message("+        Heatmaps                   ")
message("+-------------------------------------------------------------------------------")
setwd(ResDir)


genes_TGFBeta_WP <- read.csv("WPTGFbsignalling.csv")
genes_TGFBeta_GO <- read.csv("GO0007179tgfb.csv")
genes_TGFBeta_GO_sel <- read.csv("GO0017015genes.csv")
genes_TGFBeta_GO_sel <- genes_TGFBeta_GO_sel$SYMBOL

#selected_genes <- genes_TGFBeta_GO_sel
#selected_genes <- genes_TGFBeta_GO$SYMBOL
selected_genes <- genes_TGFBeta_WP$shared.name

#selected_genes <- unique(c(as.character(genes_TGFBeta_GO$SYMBOL), as.character(genes_TGFBeta_WP$shared.name), as.character(genes_TGFBeta_GO_sel )))

selected_genes <- marker_genes_0.6



Idents(matrix.su_T24) <- matrix.su_T24@meta.data$SCT_snn_res.0.4
mat_clusterMeans <- AverageExpression(matrix.su_T24, assays = "SCT", slot = "scale.data")
#mat_clusterMeans <- AverageExpression(matrix.su_T24, assays = "SCT", slot = "data")
#mat_clusterMeans[[1]] <- scale(mat_clusterMeans[[1]])

mat_clusterMeans_SCT <- mat_clusterMeans[[1]][rownames(mat_clusterMeans[[1]]) %in% selected_genes,]
mat_clusterMeans_SCT[1:5, 1:5]

mat_clusterMeans_SCT <- mat_clusterMeans_SCT[,c("0","1","2","3","4","5","9")]

mat_clusterMeans_SCT <- mat_clusterMeans_SCT - rowMeans(mat_clusterMeans_SCT)
MinVal <- min(mat_clusterMeans_SCT)
MaxVal <- max(mat_clusterMeans_SCT)



BinCols <- c( "#008837","#a6dba0","#f7f7f7","#c2a5cf","#762a83")

#MedianVal <- median(matrix_for_ht) 
#MinVal <- min(matrix_for_ht) 
#MaxVal <- max(matrix_for_ht) 

MinVal <- -2
MaxVal <- 2

steps <- ( MaxVal- MinVal) / length(BinCols)
MinVal + length(BinCols)*steps == MaxVal 

Bins <- vector()
count <- 0
for (i in seq_along(1:length(BinCols))){
  if (i == 1) {
    Bins[[i]] <-  min(mat_clusterMeans_SCT)
    count <- count + 1
  } else {
    Bins[[i]] <-  Bins[[(i-1)]] + steps
    count <- count + 1
  }
}



f3 = colorRamp2( Bins, BinCols, space = "LAB") 
f3 = colorRamp2( Bins, BinCols, space = "LAB") 

ht3 = Heatmap(t(mat_clusterMeans_SCT),  col = f3, row_title = "", column_title = "T24: Marker genes", show_row_names = TRUE, heatmap_legend_param = list(title = "MeanCentred Expr", legend_height = unit(5, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , row_names_side ="left",  show_heatmap_legend = TRUE)
ht3

#ht3_0.6 <- ht3

pdf(paste( Project, "ComplexHeatmap_clusterAverages", "TGFBeta",  "WikiPathways", "res.0.4_RmMixedCluster", ".pdf", sep="_"), width=6, height=5)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


#png(paste( Project, "DotPlot_", "TGFBeta",  "GO0017015genes.csv", "res.0.4", ".png", sep="_"), width=2000, height=500, type = "cairo")
#par(bg=NA)
#DotPlot(matrix.su_T24, features = selected_genes) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#dev.off()



#### choose option: scale.data or scale myself!!!

matrix_for_ht <- GetAssayData(matrix.su_T24, assay = "SCT", slot = "scale.data")

#matrix_for_ht <- GetAssayData(matrix.su_T24, assay = "SCT", slot = "data")
#dim(matrix_for_ht)
#matrix_for_ht <- scale(matrix_for_ht)


matrix_for_ht[1:5,1:5]
matrix_for_ht <- matrix_for_ht[rownames(matrix_for_ht) %in% selected_genes,]

split <- as.data.frame(colnames(matrix_for_ht))
colnames(split) <- "cells"
split$cluster <- matrix.su_T24@meta.data[match(split$cells, rownames(matrix.su_T24@meta.data)),]$SCT_snn_res.0.4
matrix_for_ht <- t(as.matrix(matrix_for_ht))
matrix_for_ht <- matrix_for_ht - rowMeans(matrix_for_ht)
dim(matrix_for_ht)


#00441b. Dark green
#1b7837
#5aae61
#a6dba0
#d9f0d3
#f7f7f7 whitish
#e7d4e8
#c2a5cf
#9970ab
#762a83
#40004b dark purple
#BinCols <- c("#00441b", "#1b7837", "#5aae61","#a6dba0","#d9f0d3","#f7f7f7","#e7d4e8","#c2a5cf","#9970ab","#762a83","#40004b")
BinCols <- c( "#008837","#a6dba0","#f7f7f7","#c2a5cf","#762a83")

MedianVal <- median(matrix_for_ht) 
#MinVal <- min(matrix_for_ht) + 3
#MaxVal <- max(matrix_for_ht) - 15
MinVal <- -5
MaxVal <- 10

steps <- ( MaxVal- MinVal) / length(BinCols)
MinVal + length(BinCols)*steps == MaxVal 

Bins <- vector()
count <- 0
for (i in seq_along(1:length(BinCols))){
  if (i == 1) {
  Bins[[i]] <-  min(matrix_for_ht)
  count <- count + 1
  } else {
    Bins[[i]] <-  Bins[[(i-1)]] + steps
    count <- count + 1
  }
}



f3 = colorRamp2( Bins, BinCols, space = "LAB") 

ht3 = Heatmap(matrix_for_ht,  col = f3, row_title = "", column_title = "T24: TGFBeta genes WikiPathways", show_row_names = FALSE, heatmap_legend_param = list(title = "MeanCentred Expr", legend_height = unit(5, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , split = split$cluster, row_names_side ="left",  show_heatmap_legend = TRUE)
ht3

#ht3_0.6 <- ht3

png(paste( Project, "ComplexHeatmap", "TGFBeta",  "WikiPathways", "res.0.2", "x.png", sep="_"), width=500, height=1000, type = "cairo")
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


#FeaturePlot(matrix.su_T24, features = c("Saa3", "Ccl2", "Igf2", "Ybx1", "Sparc", "Atp5b"), ncol = 4)





message("--------------------------------------------------------------------------------")
message("+        cluster 3 from T24 compare with together dataset clusters                   ")
message("+-------------------------------------------------------------------------------")

T24_clust3_cells <- metaData[metaData$SCT_snn_res.0.4 == "3",]$cell
saveRDS(T24_clust3_cells, "T24_clust3_cells.Rds")
T24_clust3_cells <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g/T24_clust3_cells.Rds")




message("--------------------------------------------------------------------------------")
message("+                   plot heatmap of genes with    GO:0009986                    ")
message("+-------------------------------------------------------------------------------")
setwd(ResDir)

#marker_genes_0.25

genes_cellSurface <- read.csv("GO-0009986_CellSurface.csv")
genes_cellSurface <- genes_cellSurface[,c(1:14)]
genes_cellSurface <- unique(genes_cellSurface$SYMBOL)

genes_cellSurface_markers0.25 <- genes_cellSurface[genes_cellSurface %in% marker_genes_0.25]
genes_cellSurface_markers0.6 <- genes_cellSurface[genes_cellSurface %in% marker_genes_0.6]



pdf(paste("DotPlot_", Project, "_cellSurfaceMarkers", "res_", reso,".pdf", sep = "_"), width=13, height=5)
par(bg=NA)
DotPlot(matrix.su_T24, features = genes_cellSurface_markers0.25 , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf(paste("FeaturePlot_", Project, "_cellSurfaceMarkers_Mif_Rala_Bst2", "res_", reso,".pdf", sep = "_"), width=5, height=13)
par(bg=NA)
FeaturePlot(matrix.su_T24, features = c("Mif", "Rala", "Bst2"), ncol = 1)
dev.off()



message("--------------------------------------------------------------------------------")
message("+        Markers all R vs all I                   ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 0.25
Idents(matrix.su_T24) <-  matrix.su_T24@meta.data$Experiment
markers_t24 <- FindMarkers(matrix.su_T24, ident.1 = "I.24", ident.2= "R.24", logfc.threshold = l2fc_cutoff, verbose = TRUE)
markers_t24 <- markers_t24[order(abs(markers_t24$avg_log2FC), decreasing = TRUE),]
head(markers_t24)
markers_t24[abs(markers_t24$avg_log2FC) > 0.6,]
write.csv(markers_t24, paste(Project, "R_vsI__markers_t24.csv", sep = "_"))

markers_t24_0.6 <- subset(markers_t24, abs(markers_t24$avg_log2FC) >= 0.6)
markers_t24_0.6 <- markers_t24_0.6[order(markers_t24_0.6$avg_log2FC),]

Idents(matrix.su_T24) <-  matrix.su_T24@meta.data[[reso]]
Idents(matrix.su_T24) <- factor(Idents(matrix.su_T24), levels = c(1,2,3,9,7,8,0,4,5,6))

pdf(paste("DotPlot_", Project, "R_vsI__markers_t24", "res_", reso,".pdf", sep = "_"), width=9, height=5)
par(bg=NA)
DotPlot(matrix.su_T24, features = rownames(markers_t24_0.6) , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


message("--------------------------------------------------------------------------------")
message("+        Markers cluster 3 (I) vs all else                   ")
message("+-------------------------------------------------------------------------------")

l2fc_cutoff <- 0.25
Idents(matrix.su_T24) <-  matrix.su_T24@meta.data[[reso]]
Idents(matrix.su_T24) <- factor(Idents(matrix.su_T24), levels = c(1,2,3,9,7,8,0,4,5,6))

markers_t24 <- FindMarkers(matrix.su_T24, ident.1 = "3", logfc.threshold = l2fc_cutoff, verbose = TRUE)
markers_t24 <- markers_t24[order(abs(markers_t24$avg_log2FC), decreasing = TRUE),]
head(markers_t24)
markers_t24[abs(markers_t24$avg_log2FC) > 0.6,]
write.csv(markers_t24, paste(Project, "clust3_vs_all__markers_t24.csv", sep = "_"))

markers_t24_0.6 <- subset(markers_t24, abs(markers_t24$avg_log2FC) >= 0.6)
markers_t24_0.6 <- markers_t24_0.6[order(markers_t24_0.6$avg_log2FC),]

markers_t24_1 <- subset(markers_t24, abs(markers_t24$avg_log2FC) >= 1)
markers_t24_1 <- markers_t24_1[order(markers_t24_1$avg_log2FC),]


pdf(paste("DotPlot_", Project, "clust3_vs_all__markers_t24", "res_", reso,"l2fc0.6.pdf", sep = "_"), width=15, height=5)
par(bg=NA)
DotPlot(matrix.su_T24, features = rownames(markers_t24_0.6) , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf(paste("DotPlot_", Project, "clust3_vs_all__markers_t24", "res_", reso,"l2fc1.pdf", sep = "_"), width=5.5, height=5)
par(bg=NA)
DotPlot(matrix.su_T24, features = rownames(markers_t24_1) , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



message("--------------------------------------------------------------------------------")
message("+        SCENIC                   ")
message("+-------------------------------------------------------------------------------")











message("--------------------------------------------------------------------------------")
message("+        CellPhoneDB                   ")
message("+-------------------------------------------------------------------------------")

normData <- GetAssayData(matrix.su_T24, assay = "SCT", slot = "data")
dim(normData)
write.csv(normData, "HarmonySCT_T24_normCounts.csv")

metaData <- matrix.su_T24@meta.data[,c("SCT_snn_res.0.4", "Treatment")]
metaData$cell <- rownames(metaData)
write.csv(metaData[,c(3,1)], "HarmonySCT_T24_metaData_clust.csv")
write.csv(metaData[,c(3,2)], "HarmonySCT_T24_metaData_treatment.csv")



