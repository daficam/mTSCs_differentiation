#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()
options(bitmapType='cairo')
#httr::set_config(httr::config(ssl_verifypeer = FALSE))

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

dataset <- "HARMONY.orig_500g"
#dataset <- "HARMONY_I_500g"
#dataset <- "HARMONY_R_500g"

if (dataset == "HARMONY.orig_500g"){
  Project <- "CTR_dscj1_HARMONY_orig_500g"
  reso <- "SCT_snn_res.0.8"
  baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g"
  
  #"12" "15" "0" "18" "11" "17" "3" "2" "13" "7" "14" "8" "10" "1" "9" "5" "16" "4" "6" "19"
  setwd(baseDir)
  #matrix.su <- readRDS("HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")
  clusters <- c(0:20) # for together dataset! - cluster 20 may need to be removed !
  row.order_tempora <-c("c12",  "c15", "c00", "c03",  "c17",  "c11", "c13", "c02", "c14",    "c18",  "c07",  "c08", "c05",  "c09", "c10",  "c16", "c01", "c19", "c06" , "c04" ) # cluster 20 removed
  # for loop=T
  #row.order_monocle <-c("c12",  "c15", "c00", "c18",  "c11",  "c17", "c03",  "c02","c13",  "c07", "c14","c08", "c10","c01",  "c09", "c05",  "c16", "c04", "c06" ,   "c19" ) # cluster 20 removed
  # for loop=F
  row.order_monocle <- c( "c12" ,"c00" , "c18", "c15", "c03" , "c02" , "c11", "c17", "c09" , "c13" ,"c08" , "c19", "c05" , "c07" , "c04" , "c06" , "c16", "c10", "c14", "c01") 
  GO_matrix <- readRDS( "HARMONY.orig_500g_ScatterPlot__GO_matrix_filter_l2fc0.25MarkerGenes")
} else if (dataset == "HARMONY_I_500g"){
  Project <- "CTR_dscj1_HARMONY_I_500g"
  reso <- "SCT_snn_res.0.8"
  baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I"
  setwd(baseDir)
  #matrix.su <- readRDS("harmony_500g_matrix.su_I.Rds")
  clusters <- c(0:16) 
  clust_order <-c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0") # for I 500g res0.8
  GO_matrix <- readRDS( "HARMONY_I_500g_ScatterPlot__GO_matrix_filter_l2fc0.25MarkerGenes")
  row.order <-c("c09",  "c16", "c11", "c05",  "c01",  "c15", "c10", "c12", "c02",  "c03",  "c06",  "c07",  "c13", "c08",  "c14", "c04",  "c00"  )
} else if (dataset == "HARMONY_R_500g"){
  Project <- "CTR_dscj1_HARMONY_R_500g"
  reso <- "SCT_snn_res.0.8"
  baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I"
  setwd(baseDir)
  #matrix.su <- readRDS("harmony_500g_matrix.su_R.Rds")
  clusters <- c(0:15) 
  clust_order <-c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8") # for R 500g res0.8
  row.order <-c("c14", "c02",  "c15", "c07",  "c13", "c06",  "c03",  "c00",  "c12", "c05" , "c10", "c01",  "c04",  "c11", "c09",  "c08" )
  GO_matrix <- readRDS( "HARMONY_R_500g_ScatterPlot__GO_matrix_filter_l2fc0.25MarkerGenes")
} else { print("please specify correct dataset!!")}

database <- "WikiPathways_2019_Mouse"



message("--------------------------------------------------------------------------------")
message("+              calculate overall mean expr for each gene                        ")
message("+-------------------------------------------------------------------------------")
library(cowplot)

ASSAY <- "SCT"
SLOT <- "data"

# get averages for gene expression for each cluster
Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.0.8
mat_clusterMeans <- AverageExpression(matrix.su, assays = ASSAY, slot = SLOT)
names(mat_clusterMeans$SCT) # "0"  "1"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2"  "20" "3"  "4"  "5"  "6"  "7"  "8"  "9" 
head(mat_clusterMeans$SCT[1])

# get averages for gene expression for all dataset - a background
Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.0
mean_gene_expr <- AverageExpression(matrix.su, assays = ASSAY, slot = SLOT)
dim(mean_gene_expr[[1]])

head(mean_gene_expr[[1]])
head(mat_clusterMeans$SCT[1])
names(mat_clusterMeans$SCT) 

mean_gene_expr1 <- as.data.frame(mean_gene_expr[[1]])
mean_gene_expr1$gene <- rownames(mean_gene_expr1) 
colnames(mean_gene_expr1)[1] <- "mean_expr_all"
mean_gene_expr1$mean_expr_c00 <- as.data.frame(mat_clusterMeans$SCT["0"])[[1]]
mean_gene_expr1$mean_expr_c01 <- as.data.frame(mat_clusterMeans$SCT["1"])[[1]]
mean_gene_expr1$mean_expr_c10 <- as.data.frame(mat_clusterMeans$SCT["10"])[[1]]
mean_gene_expr1$mean_expr_c11 <- as.data.frame(mat_clusterMeans$SCT["11"])[[1]]
mean_gene_expr1$mean_expr_c12 <- as.data.frame(mat_clusterMeans$SCT["12"])[[1]]
mean_gene_expr1$mean_expr_c13 <- as.data.frame(mat_clusterMeans$SCT["13"])[[1]]
mean_gene_expr1$mean_expr_c14 <- as.data.frame(mat_clusterMeans$SCT["14"])[[1]]
mean_gene_expr1$mean_expr_c15 <- as.data.frame(mat_clusterMeans$SCT["15"])[[1]]
mean_gene_expr1$mean_expr_c16 <- as.data.frame(mat_clusterMeans$SCT["16"])[[1]]
mean_gene_expr1$mean_expr_c17 <- as.data.frame(mat_clusterMeans$SCT["17"])[[1]]
mean_gene_expr1$mean_expr_c18 <- as.data.frame(mat_clusterMeans$SCT["18"])[[1]]
mean_gene_expr1$mean_expr_c19 <- as.data.frame(mat_clusterMeans$SCT["19"])[[1]]
mean_gene_expr1$mean_expr_c02 <- as.data.frame(mat_clusterMeans$SCT["2"])[[1]]
mean_gene_expr1$mean_expr_c20 <- as.data.frame(mat_clusterMeans$SCT["20"])[[1]]
mean_gene_expr1$mean_expr_c03 <- as.data.frame(mat_clusterMeans$SCT["3"])[[1]]
mean_gene_expr1$mean_expr_c04 <- as.data.frame(mat_clusterMeans$SCT["4"])[[1]]
mean_gene_expr1$mean_expr_c05 <- as.data.frame(mat_clusterMeans$SCT["5"])[[1]]
mean_gene_expr1$mean_expr_c06 <- as.data.frame(mat_clusterMeans$SCT["6"])[[1]]
mean_gene_expr1$mean_expr_c07 <- as.data.frame(mat_clusterMeans$SCT["7"])[[1]]
mean_gene_expr1$mean_expr_c08 <- as.data.frame(mat_clusterMeans$SCT["8"])[[1]]
mean_gene_expr1$mean_expr_c09 <- as.data.frame(mat_clusterMeans$SCT["9"])[[1]]


message("--------------------------------------------------------------------------------")
message("+                               scatter plot                                    ")
message("+-------------------------------------------------------------------------------")

lm_eqn <- function(df){
  x <- df[,1]
  y <- df[,2]
  m <- lm( x ~ y, df);
  r2 = format(summary(m)$r.squared, digits = 3)
  return(r2)
}


plot_meanExpr_clust_vs_all <- function(df=mean_gene_expr1, meanExpr_all="all", cluster , ann_number = 40) {
  
  colName1 <- paste0("mean_expr_", meanExpr_all)
  colName2 <- paste0("mean_expr_", cluster)
  df_temp <- df[,colnames(df) %in% c(colName1, colName2)]
  # select cut-off for expressed genes -remove these super low expr
  df_temp <- df_temp[df_temp$mean_expr_all > 0.0001,]
  # needed to remove mean 0 expression as cannot calculate lm otherwise
  lowExpr <- ifelse(df_temp$mean_expr_all < 0.001, "verylow", "other")
  df_temp[df_temp == 0] <- 0.00001
  df_temp <- log2(df_temp)
  scatter.log2.r2    <- lm_eqn(df_temp)
  lab.log.x      <- (max(df_temp[,1]) / 10 * 5)
  lab.log.y      <- (min(df_temp[,2]) + 5)
 
  scatter.rlog_anno <- df_temp
  df_temp$ratio <- df_temp[[1]] / df_temp[[2]] 
  df_temp$anno <- "not"
  df_temp$anno <- ifelse( df_temp[[1]] > df_temp[[2]] + 1 |  df_temp[[2]] > df_temp[[1]] + 1 , "anno", "not") 
  df_temp$anno <- factor(df_temp$anno, levels = c("anno", "not"))
  df_temp$external_gene_name <- rownames(df_temp)
  df_temp$col <- ifelse(df_temp$anno == "anno", "red", "blue")
  df_temp$diff <- df_temp[,1] - df_temp[,2]
  df_temp$lowExpr <- lowExpr
  labeldata.ann <- df_temp[df_temp$anno == "anno",]
  labeldata.ann1 <- labeldata.ann[order(labeldata.ann$diff, decreasing = TRUE),]
  labeldata.ann2 <- labeldata.ann[order(labeldata.ann$diff, decreasing = FALSE),]
  labeldata.ann1 <- labeldata.ann1[labeldata.ann1$diff > 1, ]
  labeldata.ann2 <- labeldata.ann2[labeldata.ann2$diff <= -1, ]
  
  labeldata.ann1 <- labeldata.ann1[labeldata.ann1$lowExpr == "other", ]
  labeldata.ann2 <- labeldata.ann2[labeldata.ann2$lowExpr == "other", ]
  
  labeldata.ann1 <- labeldata.ann1[1:(ann_number/2), ]
  labeldata.ann2 <- labeldata.ann2[1:(ann_number/2), ]
  labeldata.ann <- rbind(labeldata.ann1, labeldata.ann2)
  
  write.csv(df_temp[df_temp$anno == "anno" & df_temp$lowExpr == "other", c(1,2,3,7,8)], file = paste0( "DF_mean_" , cluster , "_vs_MeanExprAll_", Project , ".csv"), row.names = TRUE)
  print(nrow(df_temp[df_temp$anno == "anno" & df_temp$lowExpr == "other", c(1,2,3,7,8)]) )

  message(paste0(colName1, " ::: ", colName2, " ::: ", scatter.log2.r2, " ::: ", nrow(labeldata.ann), " labeled genes"))
  elementTextSize <- 12
  
  plt.corr <- ggplot(df_temp, aes(x = (df_temp[,1]), y=(df_temp[,2]))) + 
    geom_point(size = 2, alpha=0.3 , colour = df_temp$col) + scale_y_continuous() +
    annotate("text", x=lab.log.x, y=lab.log.y, label=paste("r^2 = ", scatter.log2.r2, sep=""), colour='navyblue', size=4 ) + 
    ggtitle(paste0("Corr.plot :::  ", cluster , " vs ", meanExpr_all )) +
    theme(aspect.ratio=1) +     
    #scale_color_discrete(col= c("red","blue"))
    geom_abline(intercept = 1, slope = 1, colour = "grey") +
    geom_abline(intercept = -1, slope = 1, colour = "grey") +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    theme_classic() + 
    geom_vline(xintercept= -9.96, linetype="dashed", color = "grey", size=1) +
    geom_hline(yintercept= -9.96, linetype="dashed", color = "grey", size=1) +
    labs(y=colName2, x = colName1) + 
    geom_label_repel(data=labeldata.ann, aes(x=labeldata.ann[,1], y=labeldata.ann[,2], label=external_gene_name),colour='purple',  size=3 )  
  return(list(plt.corr, df_temp[df_temp$anno == "anno", c(1,2,3,7,8)], labeldata.ann)) # 
}


plt.c00 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c00")
plt.c01 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c01")
plt.c02 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c02")
plt.c03 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c03")
plt.c04 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c04")
plt.c05 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c05")
plt.c06 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c06")
plt.c07 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c07")
plt.c08 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c08")
plt.c09 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c09")
plt.c10 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c10")
plt.c11 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c11")
plt.c12 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c12")
plt.c13 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c13")
plt.c14 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c14")
plt.c15 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c15")
plt.c16 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c16")
plt.c17 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c17")
plt.c18 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c18")
plt.c19 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c19")
plt.c20 <- plot_meanExpr_clust_vs_all(mean_gene_expr1, "all", "c20")

# for Together
scatterplot_dfs <- list(plt.c00[[2]], plt.c01[[2]],plt.c02[[2]],plt.c03[[2]],plt.c04[[2]],plt.c05[[2]],plt.c06[[2]],plt.c07[[2]],plt.c08[[2]],plt.c09[[2]],plt.c10[[2]],plt.c11[[2]],plt.c12[[2]],plt.c13[[2]],plt.c14[[2]],plt.c15[[2]],plt.c16[[2]],plt.c17[[2]],plt.c18[[2]],plt.c19[[2]])
# for I
scatterplot_dfs <- list(plt.c00[[2]], plt.c01[[2]],plt.c02[[2]],plt.c03[[2]],plt.c04[[2]],plt.c05[[2]],plt.c06[[2]],plt.c07[[2]],plt.c08[[2]],plt.c09[[2]],plt.c10[[2]],plt.c11[[2]],plt.c12[[2]],plt.c13[[2]],plt.c14[[2]],plt.c15[[2]],plt.c16[[2]])
# for R
scatterplot_dfs <- list(plt.c00[[2]], plt.c01[[2]],plt.c02[[2]],plt.c03[[2]],plt.c04[[2]],plt.c05[[2]],plt.c06[[2]],plt.c07[[2]],plt.c08[[2]],plt.c09[[2]],plt.c10[[2]],plt.c11[[2]],plt.c12[[2]],plt.c13[[2]],plt.c14[[2]],plt.c15[[2]])


# for Together
pdf(paste(Project, "scatterplots", "cluster_vs_all_expr", "40LabGenes_vline.pdf", sep="_"), onefile=FALSE, width=20, height=25) 
par(bg=NA)
plot_grid(plt.c00[[1]], plt.c01[[1]], plt.c02[[1]], plt.c03[[1]],
          plt.c04[[1]], plt.c05[[1]], plt.c06[[1]], plt.c07[[1]], 
          plt.c08[[1]], plt.c09[[1]], plt.c10[[1]], plt.c11[[1]],
          plt.c12[[1]], plt.c13[[1]], plt.c14[[1]], plt.c15[[1]], 
          plt.c16[[1]], plt.c17[[1]], plt.c18[[1]], plt.c19[[1]], ncol = 4, nrow = 5)
dev.off()

# for I:
pdf(paste(Project, "scatterplots", "cluster_vs_all_expr", "40LabGenes_vline.pdf", sep="_"), onefile=FALSE, width=20, height=25) 
par(bg=NA)
plot_grid(plt.c00[[1]], plt.c01[[1]], plt.c02[[1]], plt.c03[[1]],
          plt.c04[[1]], plt.c05[[1]], plt.c06[[1]], plt.c07[[1]], 
          plt.c08[[1]], plt.c09[[1]], plt.c10[[1]], plt.c11[[1]],
          plt.c12[[1]], plt.c13[[1]], plt.c14[[1]], plt.c15[[1]], 
          plt.c16[[1]], ncol = 4, nrow = 5)
dev.off()

# for R:
pdf(paste(Project, "scatterplots", "cluster_vs_all_expr", "40LabGenes_vline.pdf", sep="_"), onefile=FALSE, width=20, height=25) 
par(bg=NA)
plot_grid(plt.c00[[1]], plt.c01[[1]], plt.c02[[1]], plt.c03[[1]],
          plt.c04[[1]], plt.c05[[1]], plt.c06[[1]], plt.c07[[1]], 
          plt.c08[[1]], plt.c09[[1]], plt.c10[[1]], plt.c11[[1]],
          plt.c12[[1]], plt.c13[[1]], plt.c14[[1]], plt.c15[[1]], ncol = 4, nrow = 5)
dev.off()



message("--------------------------------------------------------------------------------")
message("+                    scatterplot genes - GO analysis                            ")
message("+-------------------------------------------------------------------------------")


if (dataset == "HARMONY.orig_500g"){
  result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/Harmony.orig_500g_Markers/result_list_harmony.orig_500g_res0.8.Rds")
} else if (dataset == "HARMONY_I_500g"){
  result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_I_500g_markers/Harmony_I_only_500g_res.0.8__result_list.Rds") 
  scatterplot_dfs <- readRDS("HARMONY_I_500g_ScatterPlot__list_cluster_pathways_filter_l2fc0.25MarkerGenes")
  } else if (dataset == "HARMONY_R_500g"){
    result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_R_500g_markers/Harmony_R_only_500g_res.0.8__result_list.Rds") 
    scatterplot_dfs <- readRDS("HARMONY_R_500g_ScatterPlot__list_cluster_pathways_filter_l2fc0.25MarkerGenes")
} else { print("please specify correct dataset!!")}




l2fc_cutoff <- 0.25
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
      #tmp_df <- unique(merge(tmp_df, ensEMBL2id, by.x = "row.names", by.y = "external_gene_name", all.x = TRUE))
      tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
      tmp_df <- subset(tmp_df, abs(tmp_df$avg_logFC) > l2fc_cutoff)
      print(rownames(tmp_df)) 
      return(rownames(tmp_df))
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

marker_genes_l2fc0.25 <- unique(unlist(unlist(marker_gene_list)))
#marker_genes_l2fc0.25 <- marker_genes_l2fc0.25[-1]


head(plt.c00[[2]])

extract_UP_and_DOWN_genes_from_scatter <- function(df, marker_gene_list = marker_genes_l2fc0.25){
  df <- df[rownames(df) %in% marker_gene_list,]
  df <- df[df$lowExpr == "other",]
  UP_genes <- rownames(df[df$mean_expr_all < df[,2],])
  DOWN_genes <- rownames(df[df$mean_expr_all > df[,2],])
  return(list(UP_genes, DOWN_genes))
}

# testing:
genes_c00 <- extract_UP_and_DOWN_genes_from_scatter(plt.c00[[2]], marker_genes_l2fc0.25)
genes_c05 <- extract_UP_and_DOWN_genes_from_scatter(plt.c05[[2]], marker_genes_l2fc0.25)
genes_c08 <- extract_UP_and_DOWN_genes_from_scatter(plt.c08[[2]], marker_genes_l2fc0.25)
genes_c07 <- extract_UP_and_DOWN_genes_from_scatter(plt.c07[[2]], marker_genes_l2fc0.25)




library("enrichR")
enrichR_DB <- as.data.frame(listEnrichrDbs())
database <- "WikiPathways_2019_Mouse"
# WikiPathways_2019_Mouse, GO_Biological_Process_2018, GO_Cellular_Component_2018, GO_Molecular_Function_2018,  BioCarta_2016, KEGG_2019_Mouse

enrichR_analysis_fun <- function( UP_DOWN_list, db){
  if (length(UP_DOWN_list[[1]]) > 0){
    enrichR_RES_UP <- as.data.frame(enrichr(UP_DOWN_list[[1]], databases = db))
    enrichR_RES_UP <- subset(enrichR_RES_UP, enrichR_RES_UP[,4] < 0.05)
    }
  if (nrow(enrichR_RES_UP) > 0){
  enrichR_RES_UP$direction <- "Upregulated_genes"
  }
  if (length(UP_DOWN_list[[2]]) > 0){
    enrichR_RES_DOWN <- as.data.frame(enrichr(UP_DOWN_list[[2]], databases = db))
    enrichR_RES_DOWN <- subset(enrichR_RES_DOWN, enrichR_RES_DOWN[,4] < 0.05) 
    }
  if (nrow(enrichR_RES_DOWN) > 0){
    enrichR_RES_DOWN$direction <- "Downregulated_genes"
  }
  if (nrow(enrichR_RES_UP) > 0 & nrow(enrichR_RES_DOWN) > 0){
  enrichR_RES <- rbind(enrichR_RES_UP, enrichR_RES_DOWN)
  } else if (nrow(enrichR_RES_UP) > 0) {
    enrichR_RES <- enrichR_RES_UP
  } else if (nrow(enrichR_RES_DOWN) > 0) { 
    enrichR_RES <- enrichR_RES_DOWN
  } else {
    enrichR_RES <- data.frame()
  }
  return(enrichR_RES)
}

# testing:
pathways_c05 <- enrichR_analysis_fun(genes_c05, db = database)
pathways_c07 <- enrichR_analysis_fun(genes_c07, db = database)


# for all of the clusters:::

list_cluster_pathways <- list()
count <- 0
ListElementNames <- list()

for (i in scatterplot_dfs){
  if (count < 10){
    ListElementName <- paste0("c0", count)
  } else {
    ListElementName <- paste0("c", count) 
    }
  count <- count + 1
  print(count)
  db.term <- paste0(database, ".Term")
  genes_cluster_x <- extract_UP_and_DOWN_genes_from_scatter(i, marker_genes_l2fc0.25 )
  pathways_cluster_x <- enrichR_analysis_fun(genes_cluster_x, db = database)
  print(nrow(pathways_cluster_x))
  list_cluster_pathways[[count]] <- pathways_cluster_x
  ListElementNames[[count]] <- ListElementName
}

names(list_cluster_pathways) <- ListElementNames
#saveRDS(list_cluster_pathways, paste( dataset, "ScatterPlot__list_cluster_pathways", "filter_l2fc0.25MarkerGenes", sep = "_"))  
  

for (i in seq_along(1: length(list_cluster_pathways))){
  db.term <- paste0(database, ".Term")
  tmp_df <- list_cluster_pathways[[i]]
  cluster_name <- names(list_cluster_pathways)[[i]]
  GO_matrix_to_add <- subset(tmp_df, tmp_df[,4] < 0.05)
  GO_matrix_to_add <- GO_matrix_to_add[,c(1,4,10)]
  GO_matrix_to_add$score <- 1 - log10(GO_matrix_to_add[,2])
  GO_matrix_to_add$score <- ifelse(GO_matrix_to_add$direction == "Upregulated_genes", GO_matrix_to_add$score, GO_matrix_to_add$score* -1)
  GO_matrix_to_add <- GO_matrix_to_add[,c(1,4)]
  colnames(GO_matrix_to_add)[2] <- cluster_name
  if (i == 1){
    GO_matrix <- GO_matrix_to_add
    print(i)
  } else {
    GO_matrix <- unique(merge(GO_matrix, GO_matrix_to_add, by.x = db.term, by.y =db.term, all.x = TRUE, all.y = TRUE ))
    print(i)
  }
}

#saveRDS(GO_matrix, paste( dataset, "ScatterPlot__GO_matrix", "filter_l2fc0.25MarkerGenes", sep = "_"))  



# LOAD IN MATRIX FOR PLOTTING!!!!!!!
#GO_matrix <- readRDS( "HARMONY_R_500g_ScatterPlot__GO_matrix_filter_l2fc0.25MarkerGenes")
#GO_matrix <- readRDS( "HARMONY_I_500g_ScatterPlot__GO_matrix_filter_l2fc0.25MarkerGenes")
#GO_matrix <- readRDS( "HARMONY.orig_500g_ScatterPlot__GO_matrix_filter_l2fc0.25MarkerGenes")



GO_matrix2 <- unique(GO_matrix)

# together:::
#which(grepl(-2.547930, GO_matrix2$c13))
#GO_matrix2 <- GO_matrix2[-38,]

# R:::
which(grepl(-2.570960, GO_matrix2$c12))
GO_matrix2 <- GO_matrix2[-10,]
which(grepl(-3.298228, GO_matrix2$c12))
which(grepl("G1 to S cell cycle control WP413", GO_matrix2$WikiPathways_2019_Mouse.Term))
GO_matrix2 <- GO_matrix2[-14,]
which(grepl("mRNA processing WP310", GO_matrix2$WikiPathways_2019_Mouse.Term))
#GO_matrix2 <- GO_matrix2[-24,]
which(grepl("WP3654", GO_matrix2$WikiPathways_2019_Mouse.Term))
GO_matrix2 <- GO_matrix2[-26,]
which(grepl("WP2902", GO_matrix2$WikiPathways_2019_Mouse.Term))
GO_matrix2 <- GO_matrix2[-c(32,33,34),]

# together...
GO_matrix2$WikiPathways_2019_Mouse.Term[ GO_matrix2$WikiPathways_2019_Mouse.Term == "p53 signaling WP2902"][2] <- "p53 signaling WP2902_x"
GO_matrix2 <- GO_matrix2[!GO_matrix2$WikiPathways_2019_Mouse.Term == "p53 signaling WP2902_x",]

# R...
GO_matrix2$WikiPathways_2019_Mouse.Term[ GO_matrix2$WikiPathways_2019_Mouse.Term == "mRNA processing WP310"][2] <- "mRNA processing WP310_x"


rownames(GO_matrix2) <- GO_matrix2[,1]
GO_matrix2 <- GO_matrix2[,-1]
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
rownames(GO_matrix2) <- gsub( " WP.*" ,  "" , rownames(GO_matrix2))

dim(GO_matrix2)

MinVal <- min(GO_matrix2, na.rm = TRUE)
MaxVal <- max(GO_matrix2, na.rm = TRUE)

GO_matrix2[is.na(GO_matrix2)] <- 0

GO_res_name <- paste(Project, "scatterplot_genes", "filt_l2fc0.25", database, sep = "_")

o1 = seriate(dist(as.matrix(GO_matrix2)), method = "GW")
o2 = seriate(dist(t(as.matrix(GO_matrix2))), method = "GW")


selectedTerms <-c("Electron Transport Chain",  "mRNA processing",  "Cytoplasmic Ribosomal Proteins", "Oxidative phosphorylation" ) #,  "mRNA processing WP310_x" ,"p53 signaling WP2902_x" ,"PluriNetWork WP1763"  , "p53 signaling WP2902" "EGFR1 Signaling Pathway WP572",

# R:::
#row.order <-c("c14", "c02",  "c15", "c07",  "c13", "c06",  "c03",  "c00",  "c12", "c05" , "c10", "c01",  "c04",  "c11", "c09",  "c08" )
# I:::
#row.order <-c("c09",  "c16", "c11", "c05",  "c01",  "c15", "c10", "c12", "c02",  "c03",  "c06",  "c07",  "c13", "c08",  "c14", "c04",  "c00"  )
# together orig:::
#row.order <-c("c12",  "c15", "c00", "c03",  "c17",  "c11", "c13", "c02", "c14",    "c18",  "c07",  "c08", "c05",  "c09", "c10",  "c16", "c01", "c19", "c06" , "c04" ) # cluster 20 removed

#row.order <- row.order_tempora
row.order <- row.order_monocle

GO_matrix3 <- GO_matrix2[ , match(row.order, colnames(GO_matrix2))]

f1 = colorRamp2( c(MinVal, MinVal/4, 0, MaxVal/4, MaxVal), c("blue4", "blue1", "white", "hotpink2", "hotpink4"), space = "RGB") 

ht1 = Heatmap(as.matrix(GO_matrix3),  col = f1, name = "",  row_title = "", column_title = "", show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , row_names_side ="left", width = unit(ncol(GO_matrix3), "cm"),row_dend_reorder = FALSE, row_dend_side = "right") # 
ht1


pdf(paste(  "ComplexHeatmap", GO_res_name, "AllTerms_MonocleOrder.pdf", sep="_"), onefile=FALSE, width=15, height=12 )
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


GO_matrix3 <- GO_matrix3[rownames(GO_matrix3) %in%  selectedTerms,]
GO_matrix3 <- GO_matrix3[match(selectedTerms, rownames(GO_matrix3)), ]

ht1 = Heatmap(as.matrix(GO_matrix3),  col = f1, name = "",  row_title = "", column_title = "", show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , row_names_side ="left", width = unit(ncol(GO_matrix3), "cm"),row_dend_reorder = FALSE, row_dend_side = "right") # 
ht1

pdf(paste(  "ComplexHeatmap", GO_res_name, "selectedTerms_TemporaOrder.pdf", sep="_"), onefile=FALSE, width=15, height=2.5 )
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()





ht3 = Heatmap(as.matrix(GO_matrix3),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", width = unit(ncol(GO_matrix3), "cm"),row_dend_reorder = FALSE, row_dend_side = "right") # 
ht3


# together: width=18, height=20 
# I: width=18, height=15 
# R: width=18, height=12 

pdf(paste(  "ComplexHeatmap", GO_res_name, "selectedTerms_nonClust.pdf", sep="_"), onefile=FALSE, width=12, height=3 )
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()





scatterplot_dfs


# with seriation::
ht2 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), row_names_side ="left", width = unit(ncol(GO_matrix2), "cm"), cluster_rows = as.dendrogram(o1[[1]]),  cluster_columns = as.dendrogram(o2[[1]]), row_dend_side = "right")
ht2

# together: width=18, height=20 
# I: width=18, height=15 
# R: width=18, height=12 








message("--------------------------------------------------------------------------------")
message("+              ecdf like plot - retired idea                        ")
message("+-------------------------------------------------------------------------------")

mean_gene_expr_ord <- as.data.frame(mean_gene_expr1[order(mean_gene_expr1$mean_expr_all, decreasing = FALSE),])
mean_gene_expr_ord$ordered <-  seq(from = 0, to = 1, length.out =  nrow(mean_gene_expr_ord))
# set cutoff for expressed genes ::: 
mean_gene_expr_ord_onlyExpr <- mean_gene_expr_ord[mean_gene_expr_ord$mean_expr_all > 0.1,] 

ggplot(mean_gene_expr_ord_onlyExpr, aes(x=log(mean_expr_all), y=ordered)) + geom_point()

# melt to show all clusters on plot:
allMeanExprData_molten <- reshape2::melt(mean_gene_expr_ord_onlyExpr, id.vars = c("gene","ordered"), value.name = "meanExpression")

unique(allMeanExprData_molten$variable)

drawECDFplot_for_1_cluster <- function(allMeanExprData_molten, cluster) {
  # cluster name in format: "c00" /  "c15"
  selected_cluster <- paste0("mean_expr_", cluster)
  df <- allMeanExprData_molten[allMeanExprData_molten$variable %in% c("mean_expr_all", selected_cluster),]
  df_cast <- reshape2::dcast(df, formula = gene ~ variable)
  df_cast$distFromMean <- df_cast$mean_expr_c15 - df_cast$mean_expr_all
  df_cast <- df_cast[abs(df_cast$distFromMean) > 0.3,]
  df_cast$metric <- abs(df_cast$distFromMean)/df_cast$mean_expr_all
  df_cast$metric <- ifelse( df_cast$metric < 1, 1/df_cast$metric, df_cast$metric)
  max(df_cast$mean_expr_all) # 11.5
  df_cast$mark_gene <- ifelse(df_cast$distFromMean > 2*df_cast$mean_expr_all | abs(df_cast$distFromMean) > 1/2*df_cast$mean_expr_all, "TRUE", "FALSE")
  print(table(df_cast$mark_gene))
  df_cast <- df_cast[order(abs(df_cast$distFromMean), decreasing = TRUE),]
  anno_genes <-  df_cast[df_cast$mark_gene == TRUE,]$gene
  #anno_genes <-  df_cast$gene[1:100]
  #anno_genes <- df_cast[abs(df_cast$distFromMean) > 1 & df_cast$metric > 1.4,]$gene
  ggplot(df, aes(x=log(meanExpression), y=ordered, color=variable)) + geom_point() + scale_color_manual(values=c('#E69F00', '#56B4E9')) + scale_size_manual(values=c(2,0.5))  + geom_text_repel(data=subset(df, gene %in% anno_genes & variable != "mean_expr_all"), aes(label=gene), show.legend=F, force=1, size=3, colour="black")
}

drawECDFplot_for_1_cluster(allMeanExprData_molten, "c15")




dev.off()


message("--------------------------------------------------------------------------------")
message("+              calculate mean expr for each gene per cluster                    ")
message("+-------------------------------------------------------------------------------")




message("--------------------------------------------------------------------------------")
message("+        investigate genes behind both direction pathway change...              ")
message("+-------------------------------------------------------------------------------")



for (i in 1:length(scatterplot_dfs)){
  scatterplot_dfs[[i]] <- scatterplot_dfs[[i]][ scatterplot_dfs[[i]]$WikiPathways_2019_Mouse.Term %in% selectedTerms, ]
}

tmp_df <- scatterplot_dfs[[16]]

mrna_process_up <- tmp_df[tmp_df$WikiPathways_2019_Mouse.Term == "mRNA processing WP310" & tmp_df$direction == "Upregulated_genes",]$WikiPathways_2019_Mouse.Genes
mrna_process_up <- unlist(strsplit(mrna_process_up, ";"))

mrna_process_down <- tmp_df[tmp_df$WikiPathways_2019_Mouse.Term == "mRNA processing WP310" & tmp_df$direction == "Downregulated_genes",]$WikiPathways_2019_Mouse.Genes
mrna_process_down <- unlist(strsplit(mrna_process_down, ";"))



head(matrix.su@meta.data,2)
table(matrix.su@meta.data[,c("Experiment", "SCT_snn_res.0.8")])




