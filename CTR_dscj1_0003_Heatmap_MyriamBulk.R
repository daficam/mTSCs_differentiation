#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()


library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Seurat")  # seurat v3.02
library("dplyr")
theme_set(theme_cowplot())
library(circlize)
library(ComplexHeatmap)
library(matrixStats)


baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I"
setwd(baseDir)

sampleTable <- read.csv("sampleTable_MERGED_53.csv")


dataset <- "I"
#dataset <- "R"




#load matrix of R or I
if(dataset == "I"){
  Project <- "CTR_dscj1_harmony_I_500g_MyriamBulk"
  matrix.su <- readRDS("HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_I.Rds")
  dataset_name <- "Inhibit_32_500g"
  #matrix.umap <- readRDS("Harmony_Inhibit_32_500g_matrix.umap_.Rds")
  cluster_order_mean <- c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0")
} else if(dataset == "R"){
  Project <- "CTR_dscj1_harmony_R_500g_MyriamBulk"
  matrix.su <- readRDS("HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_R.Rds")
  dataset_name <- "Remove_29_500g"
  #matrix.umap <- readRDS("Harmony_Remove_29_500g_matrix.umap_.Rds")
  cluster_order_mean <- c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8")
} else {
  print("no dataset assigned")
}
# order of clusters by time: add sth here for order



bulk_troph  <- read.csv("MyriamTrophBulk_NatComm_BulkRNA_normCounts2.csv") 
#-> all secreted proteins in the bulk trophoblast + expression data (maybe start with this list)
#  cut -d "," -f1-39 MyriamTrophBulk_NatComm_BulkRNA_normCounts.csv > MyriamTrophBulk_NatComm_BulkRNA_normCounts2.csv

colnames(bulk_troph)
bulk_troph <- bulk_troph[,-c(2:15,24:36)]


message("--------------------------------------------------------------------------------")
message("+                   Dot plot for chosen markers BY CLUSTER                      ")
message("+-------------------------------------------------------------------------------")


# show gene expression per cluster 
Idents(matrix.su) <-  matrix.su@meta.data$SCT_snn_res.0.8
Idents(matrix.su) <- factor(Idents(matrix.su), levels = cluster_order_mean )


#divide genes into 3 groups- high in early, high in late, intermediate.

#do heatmap to see if cluster specific secretome
#come back to doplot


# -------------------   get matrix.su expression matrix --- SCT  data here ----

expr_mat_per_cluster <- AverageExpression(matrix.su, assays = "SCT", slot = "scale.data")


rv <- rowVars(as.matrix(expr_mat_per_cluster[[1]]))
select <- order(rv, decreasing = TRUE)[seq_len(min(3000, length(rv)))]

expr_mat_per_cluster <- scale(expr_mat_per_cluster[[1]])
expr_mat_per_cluster <- expr_mat_per_cluster[select,]









bulk_troph <- bulk_troph[!duplicated(bulk_troph$external_gene_name),]
rownames(bulk_troph) <- bulk_troph$external_gene_name
#mat_bulk <- secr_troph[, c(3:16)]
#mat_bulk <- secr_troph[, c(3:36)]
mat_bulk <- bulk_troph[, c(2:9)]
mat_bulk_scaled <- as.data.frame(t(scale(t(mat_bulk))))
mat_bulk_scaled <- na.omit(mat_bulk_scaled)

min(mat_bulk)
max(mat_bulk)
min(mat_bulk_scaled, na.rm = T)
max(mat_bulk_scaled, na.rm = T)
min(expr_mat_per_cluster, na.rm = T)
max(expr_mat_per_cluster, na.rm = T)



#mat_scaled <- means_tab2[rownames(means_tab2) %in% rownames(mat_bulk_scaled),]
mat_scaled <- expr_mat_per_cluster[rownames(expr_mat_per_cluster) %in% rownames(mat_bulk_scaled),]


mat_bulk_scaled <- mat_bulk_scaled[match( rownames(mat_scaled), rownames(mat_bulk_scaled)),]
rownames(mat_bulk_scaled) == rownames(mat_scaled)





MeansExpr <- expr_mat_per_cluster


### for MEAN expresion data::: 
#rownames(MeansExpr) <- MeansExpr$gene
#MeansExpr <- MeansExpr[order(MeansExpr$meanExpr, decreasing = TRUE),]
MeansExpr_top10pc <- MeansExpr[c(1:(nrow(MeansExpr)/10)),]
MeansExpr_top100 <- MeansExpr[c(1:100),]
MeansExpr_top50 <- MeansExpr[c(1:50),]

mat_scaled2 <- mat_scaled[rownames(mat_scaled) %in% rownames(MeansExpr_top100) , ]
mat_bulk_scaled2 <- mat_bulk_scaled[match( rownames(mat_scaled2), rownames(mat_bulk_scaled)),]
rownames(mat_bulk_scaled2) == rownames(mat_scaled2)

head(mat_bulk_scaled2,3)
#PL_024
#PL_X3
#PL_X24
#in this exact order (STEM->Diff).
mat_bulk_scaled2 <- mat_bulk_scaled2[,c( "lane5_PL_024_1_L005", "lane5_PL_024_2_L005", "lane5_PL_X3_1_L005",  "lane5_PL_X3_2_L005" ,   "lane5_PL_X24_1_L005", "lane5_PL_X24_2_L005")]
colnames(mat_bulk_scaled2) <- c("TSC_1", "TSC_2", "TSC_3h_MEKi_1", "TSC_3h_MEKi_2", "TSC_24h_MEKi_1", "TSC_24h_MEKi_2")

f1 = colorRamp2( c(min(mat_scaled2, na.rm = T),0, max(mat_scaled2, na.rm = T)/2), c("green3", "grey95",  "darkorchid4"), space = "LAB") # "darkgreen", "green3", "grey95", "slateblue", "darkorchid4"
f2 = colorRamp2( c(min(mat_bulk_scaled2, na.rm = T ) ,0, max(mat_bulk_scaled2, na.rm = T)), c("green3",  "grey95", "darkorchid4"), space = "LAB")
lgd1 = Legend(col_fun = f1, title = "scRNA-seq", at = c(min(mat_scaled2, na.rm = T),  0, max(mat_scaled2, na.rm = T))  )
lgd2 = Legend(col_fun = f2, title = "bulk RNA-seq",  at = c(min(mat_bulk_scaled2, na.rm = T), 0, max(mat_bulk_scaled2, na.rm = T)) )
#pd = packLegend(lgd1, lgd2, direction = "vertical", max_width = unit(3, "cm"))

ht1 = Heatmap(as.matrix(mat_scaled2),  col = f1, name = "sc RNA-seq",  row_title = "", column_title = "sc RNA-seq", show_row_names = TRUE, heatmap_legend_param = list(title = "scRNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(10, "cm"), row_names_side ="left")

ht2 = Heatmap(as.matrix(mat_bulk_scaled2),  col = f2, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "bulk RNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, width = unit(6, "cm"), row_names_side ="left")

ht_list = ht1 + ht2
ht_list

pdf(paste(dataset, "_Fig_", "Htmap_", "_bulk_troph__top100MVGs", "6placSamples", "SCT_slot_scale.data_scaled.pdf", sep=""), onefile=FALSE, width=15, height=20) # 12 for 50 mvg, 20 for 100mvg 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()





### for MAX MEANS::: 
MeansExpr <- Max_Means
rownames(MeansExpr) <- MeansExpr$gene
MeansExpr_sc_secrTroph <- as.data.frame(MeansExpr[MeansExpr$gene %in%  rownames(secr_troph),])
MeansExpr_sc_secrTroph <- MeansExpr_sc_secrTroph[order(MeansExpr_sc_secrTroph$max_means, decreasing = TRUE),]
MeansExpr_sc_secrTroph_top10pc <- MeansExpr_sc_secrTroph[c(1:(nrow(MeansExpr_sc_secrTroph)/10)),]

mat_scaled2 <- mat_scaled[rownames(mat_scaled) %in% MeansExpr_sc_secrTroph_top10pc$gene , ]
mat_bulk_scaled2 <- mat_bulk_scaled[match( rownames(mat_scaled2), rownames(mat_bulk_scaled)),]
rownames(mat_bulk_scaled2) == rownames(mat_scaled2)

f1 = colorRamp2( c(min(mat_scaled2, na.rm = T),0, max(mat_scaled2, na.rm = T)), c("green3", "grey95",  "darkorchid4"), space = "LAB") # "darkgreen", "green3", "grey95", "slateblue", "darkorchid4"
f2 = colorRamp2( c(min(mat_bulk_scaled2, na.rm = T ) ,0,max(mat_bulk_scaled2, na.rm = T)), c("green3",  "grey95", "darkorchid4"), space = "LAB")
lgd1 = Legend(col_fun = f1, title = "scRNA-seq", at = c(min(mat_scaled2, na.rm = T),  0, max(mat_scaled2, na.rm = T))  )
lgd2 = Legend(col_fun = f2, title = "bulk RNA-seq",  at = c(min(mat_bulk_scaled2, na.rm = T), 0, max(mat_bulk_scaled2, na.rm = T)) )
pd = packLegend(lgd1, lgd2, direction = "vertical", max_width = unit(3, "cm"))

ht1 = Heatmap(as.matrix(mat_scaled2),  col = f1, name = "sc RNA-seq",  row_title = "", column_title = "sc RNA-seq", show_row_names = TRUE, heatmap_legend_param = list(title = "scRNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(10, "cm"), row_names_side ="left")
ht2 = Heatmap(as.matrix(mat_bulk_scaled2),  col = f2, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "bulk RNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, width = unit(6, "cm"), row_names_side ="left")

ht_list = ht1 + ht2
ht_list

pdf(paste(dataset, "_Fig_", "Htmap_", "_bulk_NEWsecretome_troph__top10pc_in_sc_", "8placSamples", "Integrated_slot_scale.data_scaled_MAX_MEANS.pdf", sep=""), onefile=FALSE, width=15, height=13.5) 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()
























# -------------------   get matrix.su expression matrix --- integrated scale.data here ----

#expr_mat <- GetAssayData(matrix.su, assay = "RNA", slot = "data")
#expr_mat <- GetAssayData(matrix.su, assay = "integrated", slot = "scale.data")
#all_feats <- rownames(expr_mat)

DefaultAssay(matrix.su) <- "integrated"
#data.featsx <- FetchData(object = matrix.su, vars = all_feats, slot = "data") 

#For a heatmap or dotplot of markers, the scale.data in the RNA assay should be used.???

data.feats$id <- Idents(matrix.su)
id.levels <- levels(x = data.feats$id)
data.feats$id <- as.vector(x = data.feats$id) 
data.feats[1:5,1:5]
min(data.feats[,c(1: (ncol(data.feats) - 1) )])

data.means <- apply(
  X = unique(x = data.feats$id),
  FUN = function(ident) {
    data.use <- data.feats[data.feats$id == ident, 1:(ncol(x = data.feats) - 1), drop = FALSE]
    avg.exp <- apply(
      X = data.use,
      MARGIN = 2,
      FUN = function(x) {
        #return(mean(x = expm1(x = x)))
        return(mean( x))
      }
    )
    return(avg.exp = avg.exp)
  }
)

names(x = data.means) <- unique(x = data.feats$id)
data.means <- lapply(
  X = names(x = data.means),
  FUN = function(x) {
    data.use <- as.data.frame(x = data.means[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  }
)

data.means <- do.call(what = 'rbind', args = data.means)
if (!is.null(x = id.levels)) {
  data.means$id <- factor(x = data.means$id, levels = id.levels)
}


dim(data.means)
head(data.means)
colnames(data.means) <- c("mean_expr", "gene", "cluster_id")
max(data.means[,c(1)])
min(data.means[,c(1)])

data.means[data.means$mean_expr == max(data.means[,c(1)]),] # Kcnh321, clust 23

data.means2[rownames(data.means2) == "Kcnh3"  ,] # Kcnh321, clust 23



library(tidyr)
data.means2 <- spread(data.means, cluster_id, mean_expr)
#data.means2
head(data.means2,4)
rownames(data.means2) <- data.means2$gene
data.means2 <- data.means2[,-1]
min(data.means2)
max(data.means2)

#saveRDS(data.means2, "ExprMat_I_assayIntegrated_slot_scale.Data_cluster_data.means.rds")
#saveRDS(data.means2, "ExprMat_R_assayIntegrated_slot_scale.Data_cluster_data.means.rds")
#saveRDS(data.means2, "ExprMat_R_assayIntegrated_slot_scale.Data_cluster_data.means_exp1.rds")


# -------------------   load data.means object to plot heatmap ---------------------

#data.means2 <- readRDS("ExprMat_R_assayIntegrated_slot_scale.Data_cluster_data.means_exp1.rds")
#data.means2 <- readRDS("ExprMat_R_assayIntegrated_slot_scale.Data_cluster_data.means.rds")
#data.means2 <- readRDS("ExprMat_I_assayIntegrated_slot_scale.Data_cluster_data.means_exp1.rds")
data.means2 <- readRDS("ExprMat_I_assayIntegrated_slot_scale.Data_cluster_data.means.rds")

rownames(data.means2) <- gsub( "integrated_", "", rownames(data.means2))






data.means2_scaled <- as.data.frame(t(scale(t(data.means2))))
min(data.means2_scaled, na.rm = T)
max(data.means2_scaled, na.rm = T)

min(data.means2)
max(data.means2)

data.means2_scaled <- na.omit(data.means2_scaled)

#rownames(secr_troph_top10pc) <- secr_troph_top10pc$external_gene_name
#mat_bulk <- secr_troph_top10pc[, c(3:16)]
#mat_bulk_scaled <- as.data.frame(t(scale(t(mat_bulk))))

secr_troph <-secr_troph[order(secr_troph$chromosome_name),]
secr_troph <- secr_troph[!duplicated(secr_troph$external_gene_name),]
rownames(secr_troph) <- secr_troph$external_gene_name
mat_bulk <- secr_troph[, c(3:16)]
mat_bulk_scaled <- as.data.frame(t(scale(t(mat_bulk))))
mat_bulk_scaled <- na.omit(mat_bulk_scaled)

min(mat_bulk)
max(mat_bulk)
min(mat_bulk_scaled, na.rm = T)
max(mat_bulk_scaled, na.rm = T)



mat_scaled <- data.means2[rownames(data.means2) %in% rownames(mat_bulk_scaled),]
#mat_scaled <- data.means2_scaled[rownames(data.means2_scaled) %in% rownames(mat_bulk_scaled),]

#mat_bulk <- secr_troph_top10pc[secr_troph_top10pc$external_gene_name %in% rownames(mat), c(3:16,19)]
#mat_bulk <- secr_troph[secr_troph$external_gene_name %in% rownames(mat), c(3:16,19)]
#rownames(mat_bulk) <- mat_bulk$external_gene_name
#mat_bulk <- mat_bulk[,-ncol(mat_bulk)]
#mat_bulk_scaled <- as.data.frame(t(scale(t(mat_bulk))))
#mat_scaled <- as.data.frame(t(scale(t(mat))))


mat_bulk_scaled <- mat_bulk_scaled[match( rownames(mat_scaled), rownames(mat_bulk_scaled)),]
rownames(mat_bulk_scaled) == rownames(mat_scaled)

f1 = colorRamp2( c(min(mat_scaled),0,max(mat_scaled)), c("green3", "grey95",  "darkorchid4"), space = "LAB") # "darkgreen", "green3", "grey95", "slateblue", "darkorchid4"
f2 = colorRamp2( c(min(mat_bulk_scaled, na.rm = T ) ,0,max(mat_bulk_scaled, na.rm = T)), c("green3",  "grey95", "darkorchid4"), space = "LAB")
lgd1 = Legend(col_fun = f1, title = "scRNA-seq", at = c(min(mat_scaled),  0, max(mat_scaled))  )
lgd2 = Legend(col_fun = f2, title = "bulk RNA-seq",  at = c(min(mat_bulk_scaled, na.rm = T), 0, max(mat_bulk_scaled, na.rm = T)) )
pd = packLegend(lgd1, lgd2, direction = "vertical")

ht1 = Heatmap(as.matrix(mat_scaled),  col = f1, name = "sc RNA-seq",  row_title = "", column_title = "sc RNA-seq", show_row_names = TRUE, heatmap_legend_param = list(title = "scRNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(10, "cm"), row_names_side ="left")
ht2 = Heatmap(as.matrix(mat_bulk_scaled),  col = f2, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "bulk RNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, width = unit(6, "cm"), row_names_side ="left")

ht_list = ht1 + ht2
ht_list

heatmap_height <- nrow(mat_scaled)*0.2+1
#pdf(paste("Fig_xx_", "ComplexHeatmap_", "dataset_",dataset, "_bulk_secretome_troph_all", "v2_scale.data_exp1_scaled.pdf", sep=""), onefile=FALSE, width=12, height=15) 
pdf(paste("Fig_xx_", "ComplexHeatmap_", "dataset_",dataset, "_bulk_secretome_troph_all_", "v2_scale.data_notscaled.pdf", sep=""), onefile=FALSE, width=12, height=15) 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()






###################################################################################################################
Expr_I_data <- GetAssayData(matrix.su, assay = "integrated", slot = "data")
dim(Expr_I_data)
Expr_I_data[1:5,1:5]
Means_1 <- as.data.frame(rowMeans(as.data.frame(Expr_I_data[c(1:10000),])))
Means_2 <- as.data.frame(rowMeans(as.data.frame(Expr_I_data[c(10001:nrow(Expr_I_data)),])))
colnames(Means_1) <- "meanExpr"
colnames(Means_2) <- "meanExpr"
Means_1$gene <- rownames(Means_1)
Means_2$gene <- rownames(Means_2)
MeansExpr_I_data <- rbind(Means_1, Means_2)
saveRDS(MeansExpr_I_data, "MeansExpr_I_integrated.data.rds")


Expr_R_data <- GetAssayData(matrix.su, assay = "integrated", slot = "data")
dim(Expr_R_data)
Expr_R_data[1:5,1:5]
Means_1 <- as.data.frame(rowMeans(as.data.frame(Expr_R_data[c(1:10000),])))
Means_2 <- as.data.frame(rowMeans(as.data.frame(Expr_R_data[c(10001:nrow(Expr_R_data)),])))
colnames(Means_1) <- "meanExpr"
colnames(Means_2) <- "meanExpr"
Means_1$gene <- rownames(Means_1)
Means_2$gene <- rownames(Means_2)
MeansExpr_R_data <- rbind(Means_1, Means_2)
saveRDS(MeansExpr_R_data, "MeansExpr_R_integrated.data.rds")
###################################################################################################################

MeansExpr <- readRDS( "MeansExpr_I_integrated.data.rds")
MeansExpr <- readRDS( "MeansExpr_R_integrated.data.rds")

#MeansExpr <- MeansExpr_I_data
#MeansExpr <- MeansExpr_R_data

rownames(MeansExpr) <- MeansExpr$gene

MeansExpr_sc_secrTroph <- as.data.frame(MeansExpr[MeansExpr$gene %in%  rownames(secr_troph),])
MeansExpr_sc_secrTroph <- MeansExpr_sc_secrTroph[order(MeansExpr_sc_secrTroph$meanExpr, decreasing = TRUE),]

MeansExpr_sc_secrTroph_top10pc <- MeansExpr_sc_secrTroph[c(1:(nrow(MeansExpr_sc_secrTroph)/10)),]






data.means2_scaled <- as.data.frame(t(scale(t(data.means2))))
min(data.means2_scaled, na.rm = T)
max(data.means2_scaled, na.rm = T)
min(data.means2)
max(data.means2)

data.means2_scaled <- na.omit(data.means2_scaled)

data.means2_scaled <- data.means2_scaled[rownames(data.means2_scaled) %in% MeansExpr_sc_secrTroph_top10pc$gene , ]

data.means2x <- data.means2[rownames(data.means2) %in% MeansExpr_sc_secrTroph_top10pc$gene , ]
min(data.means2x)
max(data.means2x)




secr_troph <-secr_troph[order(secr_troph$chromosome_name),]
secr_troph <- secr_troph[!duplicated(secr_troph$external_gene_name),]
rownames(secr_troph) <- secr_troph$external_gene_name
mat_bulk <- secr_troph[, c(3:16)]
mat_bulk_scaled <- as.data.frame(t(scale(t(mat_bulk))))
mat_bulk_scaled <- na.omit(mat_bulk_scaled)

min(mat_bulk)
max(mat_bulk)
min(mat_bulk_scaled, na.rm = T)
max(mat_bulk_scaled, na.rm = T)


mat_scaled <- data.means2x[rownames(data.means2x) %in% rownames(mat_bulk_scaled),]
#mat_scaled <- data.means2_scaled[rownames(data.means2_scaled) %in% rownames(mat_bulk_scaled),]

mat_bulk_scaled <- mat_bulk_scaled[match( rownames(mat_scaled), rownames(mat_bulk_scaled)),]
rownames(mat_bulk_scaled) == rownames(mat_scaled)

f1 = colorRamp2( c(min(mat_scaled, na.rm = T),0,max(mat_scaled, na.rm = T)), c("green3", "grey95",  "darkorchid4"), space = "LAB") # "darkgreen", "green3", "grey95", "slateblue", "darkorchid4"
f2 = colorRamp2( c(min(mat_bulk_scaled, na.rm = T ) ,0,max(mat_bulk_scaled, na.rm = T)), c("green3",  "grey95", "darkorchid4"), space = "LAB")
lgd1 = Legend(col_fun = f1, title = "scRNA-seq", at = c(min(mat_scaled, na.rm = T),  0, max(mat_scaled, na.rm = T))  )
lgd2 = Legend(col_fun = f2, title = "bulk RNA-seq",  at = c(min(mat_bulk_scaled, na.rm = T), 0, max(mat_bulk_scaled, na.rm = T)) )
pd = packLegend(lgd1, lgd2, direction = "vertical")

ht1 = Heatmap(as.matrix(mat_scaled),  col = f1, name = "sc RNA-seq",  row_title = "", column_title = "sc RNA-seq", show_row_names = TRUE, heatmap_legend_param = list(title = "scRNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(10, "cm"), row_names_side ="left")
ht2 = Heatmap(as.matrix(mat_bulk_scaled),  col = f2, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "bulk RNA-seq Expression (scaled)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, width = unit(6, "cm"), row_names_side ="left")

ht_list = ht1 + ht2
ht_list

#pdf(paste("Fig_xx_", Project, "ComplexHeatmap_", "dataset_",dataset, "_bulk_secretome_troph__top10pc_in_sc", "_scale.data_exp1.pdf", sep=""), onefile=FALSE, width=10, height=13.5) 
pdf(paste("Fig_xx_", Project, "ComplexHeatmap_", "dataset_",dataset, "_bulk_secretome_troph__top10pc_in_sc_", "v3_scale.data_notscaled.pdf", sep=""), onefile=FALSE, width=10, height=13.5) 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


































genes_troph_secr <- rownames(secr_troph)
genes_troph_secr <- genes_troph_secr[genes_troph_secr %in% rownames(GetAssayData(matrix.su, slot = "scale.data"))]

DefaultAssay(matrix.su)

#matrix.su <- ScaleData(matrix.su, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(matrix.su)) #, features = rownames(matrix.su_I))

png(paste("Fig_xx_", "DoHeatmap_", "dataset_", dataset, "_bulk_secretome_troph_all_only_sc", ".png", sep=""), width=1200, height=2500, type = "cairo") 
par(bg=NA)
DoHeatmap(matrix.su, features =  genes_troph_secr[1:5], assay = "integrated",  angle = 90) + NoLegend()
dev.off()


doHeatMap(nbt,genes.use = markers.use,slim.col.label = TRUE,remove.key = TRUE,cexRow=0.1)
doHeatMap(nbt,genes.use = markers.use.neuronal,slim.col.label = TRUE,remove.key = TRUE,cells.use = which.cells(nbt,c(1,2,3,8)))












