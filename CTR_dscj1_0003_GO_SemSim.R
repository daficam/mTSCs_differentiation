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
  markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/Harmony.orig_500g_Markers_GO"
  setwd(markerGODir)
  GOBP_df <- read.csv("NewAlign51_SCT_HARMONY_500g__res_0.8__GOBP_matrix_l2fc0.25_.csv", row.names = "X")
  clusters <- c(0:20) # for together dataset! - cluster 20 may need to be removed !
  } else if (dataset == "HARMONY_I_500g"){
    Project <- "CTR_dscj1_HARMONY_I_500g"
    reso <- "SCT_snn_res.0.8"
    markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_I_500g_GO"
    setwd(markerGODir)
    GOBP_df <- read.csv("NewAlign51_HARMONY_Separate__Inhibit_32_500g_res.0.8_GOBP_matrix_l2fc0.25_.csv", row.names = "X")
    clusters <- c(0:16) 
  } else if (dataset == "HARMONY_R_500g"){
    Project <- "CTR_dscj1_HARMONY_R_500g"
    reso <- "SCT_snn_res.0.8"
    markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/HARMONY_SEPARATE_R_500g_GO"
    setwd(markerGODir)
    GOBP_df <- read.csv("NewAlign51_HARMONY_Separate__Remove_29_500g_res.0.8_GOBP_matrix_l2fc0.25_.csv", row.names = "X")
    clusters <- c(0:15) 
  } else { print("please specify correct dataset!!")}
  

head(GOBP_df,2)
rownames(GOBP_df) <- GOBP_df[,1]
colnames(GOBP_df) <- gsub( "X0",  "X", colnames(GOBP_df))



message("--------------------------------------------------------------------------------")
message("+               Gene cluster semantic similarity measurement                    ")
message("+-------------------------------------------------------------------------------")
# https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html

#library(biomaRt)
#ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene_id'), mart = ensembl)          

library(GOSemSim)
mmGOBP <- godata('org.Mm.eg.db', ont="BP",  keytype = "SYMBOL")
#mmGOMF <- godata('org.Mm.eg.db', ont="MF", keytype = "SYMBOL")



get_GO_for_clust_pair <- function(GOBP_df_X, first_clust, second_clust){
  if (first_clust > second_clust){ print("First cluster has to be smaller than second!!!") } 
  else {
    #if (first_clust < 10) { 
    #  first_clust <- paste0(0, first_clust)
    #}
    sel_col <- paste0("X", first_clust, ".vs.", second_clust)
    rownames(GOBP_df_X) <- gsub( ".*\\(", "", rownames(GOBP_df_X) )
    rownames(GOBP_df_X) <- gsub( "\\)", "", rownames(GOBP_df_X) )
    idx <- which(colnames(GOBP_df_X) %in% sel_col)
    GO_df_filtered <- na.omit(GOBP_df_X[, c(1, idx)])
    return(rownames(GO_df_filtered))
  }
}

# testing for single pair of GO terms cluster pair::: 
GO_terms_test1 <- get_GO_for_clust_pair(GOBP_df, first_clust = 1, second_clust = 2)
GO_terms_test2 <- get_GO_for_clust_pair(GOBP_df, first_clust = 1, second_clust = 3)
mgoSim(GO_terms_test1, GO_terms_test2, semData=mmGOBP, measure="Wang", combine="BMA")


# set up table to fill in with SCORES:
#clusters <- c(0:20) # for together dataset! - cluster 20 may need to be removed !
#clusters <- c(0:16) # for I dataset! 
#????clusters <- c(0:16) # for R dataset! 
clusterPairs <- expand.grid(clusters,clusters)
clusterPairs <- clusterPairs[clusterPairs$Var1 != clusterPairs$Var2 , ]
clusterPairs <- clusterPairs[clusterPairs$Var1 < clusterPairs$Var2 , ]
clusterPairs$name <- paste0("X", clusterPairs$Var1, ".vs.", clusterPairs$Var2)

GO_semSim_tbl <- dplyr::as_tibble(matrix(0, ncol = length(clusterPairs$name), nrow = length(clusterPairs$name) ))
rownames(GO_semSim_tbl) <- clusterPairs$name
colnames(GO_semSim_tbl) <- clusterPairs$name
GO_semSim_tbl <- as.data.frame(GO_semSim_tbl)
rownames(GO_semSim_tbl) <- clusterPairs$name
colnames(GO_semSim_tbl) <- clusterPairs$name

calculate_mgoSim_for_2_clustPairs <- function(df, clusterPair1, clusterPair2, semData=mmGOBP){
  # clusterPair1 and clusterPair2 provided as vector or 2 values in increasing order:  e.g. c(1,2)
  GO_terms_clustPair1 <- get_GO_for_clust_pair(GOBP_df, first_clust = clusterPair1[[1]], second_clust = clusterPair1[[2]])
  GO_terms_clustPair2 <- get_GO_for_clust_pair(GOBP_df, first_clust = clusterPair2[[1]], second_clust = clusterPair2[[2]])
  score <- mgoSim(GO_terms_clustPair1, GO_terms_clustPair2, semData=mmGOBP, measure="Wang", combine="BMA")
  pair1Name <- paste0("X", clusterPair1[[1]], ".vs.", clusterPair1[[2]])
  pair2Name <- paste0("X", clusterPair2[[1]], ".vs.", clusterPair2[[2]])
  idx1 <- which(rownames(GO_semSim_tbl) == pair1Name)
  idx2 <- which(colnames(GO_semSim_tbl) == pair2Name)
  #GO_semSim_tbl[idx1, idx2] <- score
  print(paste(idx1, "," ,idx2,", score: ", score))
  return((c(idx1, idx2, score)))
}

# testing :::
res_semSim <- calculate_mgoSim_for_2_clustPairs(GOBP_df, clusterPair1= c(1,2), clusterPair2= c(1,3))
GO_semSim_tbl[res_semSim[[1]], res_semSim[[2]]] <- res_semSim[[3]]
res_semSim <- calculate_mgoSim_for_2_clustPairs(GOBP_df, clusterPair1= c(4,8), clusterPair2= c(6,19))
GO_semSim_tbl[res_semSim[[1]], res_semSim[[2]]]


for (i in 1:nrow(clusterPairs)){ # 1:10) { 
  clusterPair1 <- c(clusterPairs$Var1[[i]], clusterPairs$Var2[[i]])
  #print(clusterPair1)
  for (j in 1:nrow(clusterPairs)){ # 1:10) { 
    if (j >= i) {
    clusterPair2 <- c(clusterPairs$Var1[[j]], clusterPairs$Var2[[j]])
    res_semSim <- calculate_mgoSim_for_2_clustPairs(GOBP_df, clusterPair1= clusterPair1, clusterPair2= clusterPair2)
    GO_semSim_tbl[res_semSim[[1]], res_semSim[[2]]] <- res_semSim[[3]] }
  }
}
GO_semSim_tbl[1:10,1:10]


#saveRDS(GO_semSim_tbl, "GO_semSim_tbl_I.Rds")
#saveRDS(GO_semSim_tbl, "GO_semSim_tbl_R.Rds")
#GO_semSim_tbl <- readRDS("GO_semSim_tbl_I.Rds")
#GO_semSim_tbl <- readRDS("GO_semSim_tbl_R.Rds")
GO_semSim_tbl <- readRDS("GO_semSim_tbl.Rds")

GO_semSim_tbl[1:10,1:10]
upTri <- upper.tri(GO_semSim_tbl, diag = FALSE)

GO_semSim_tbl2 <- GO_semSim_tbl
library(gdata)
lowerTriangle(GO_semSim_tbl2, diag=FALSE, byrow=TRUE) <- upperTriangle(GO_semSim_tbl2, diag=FALSE, byrow=FALSE)
GO_semSim_tbl2[1:10,1:10]


# extract T0 cluster vs rest!!
T0_clust <- 12  # 12 for together,  9 for I, 14 for R
T0_comparisons <- colnames(GO_semSim_tbl2)
T0_comparisons <- T0_comparisons[grep(T0_clust, T0_comparisons)]

GO_semSim_tbl_T0 <- GO_semSim_tbl2[rownames(GO_semSim_tbl2) %in% T0_comparisons, colnames(GO_semSim_tbl2) %in% T0_comparisons]


dim(GO_semSim_tbl_T0)
mat <- as.matrix(GO_semSim_tbl_T0)
colnames(mat) <- gsub( "X", "", colnames(mat))
rownames(mat) <- gsub( "X", "", rownames(mat))
mat[1:10,1:10]

changeColNames <- as.data.frame(colnames(mat))
colnames(changeColNames) <- "colname1"
changeColNames$clust1 <- gsub( ".vs.*", "", changeColNames$colname1)
changeColNames$clust2 <- gsub( ".*vs.", "", changeColNames$colname1)
changeColNames$clust1 <- as.numeric(as.character(changeColNames$clust1 ))
changeColNames$clust2 <- as.numeric(as.character(changeColNames$clust2 ))
changeColNames$new_colname <- ifelse(changeColNames$clust1 < T0_clust, paste0(changeColNames$clust2, ".vs.", changeColNames$clust1) , paste0(changeColNames$colname1))

colnames(mat) <- changeColNames$new_colname
rownames(mat) <- changeColNames$new_colname


f1 = colorRamp2( c(0, 0.2, 0.4, 0.5, 0.65, 0.75,  0.9, 1), c("white", "#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8","#807dba", "#6a51a3", "#4a1486"), space = "RGB")  # purple

#f1 = colorRamp2( c(0, 1), c("white", "darkblue"), space = "RGB") 
#f1 = colorRamp2( c(0, 0.4, 0.55, 0.7, 0.9, 1), c("white", "#fee0d2", "#fcbba1","#fc9272", "#de2d26", "#99000d"), space = "RGB")  # red
ht1 = Heatmap(mat,  col = f1, name = "",  row_title = "", column_title = paste0(dataset, ": GO sem sim score" ), show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "combined BMA score", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left") # width = unit(140, "cm"),
ht1

pdf(paste( Project, "GOSemSim___ComplexHeatmap_to_T0_cluster", "purple.pdf", sep="_"), onefile=FALSE, width=12, height=10)
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

ht1b = Heatmap(mat,  col = f1, name = "",  row_title = "", column_title = paste0(dataset, ": GO sem sim score" ), show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "combined BMA score", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , row_names_side ="left") # width = unit(140, "cm"),
ht1b


mat2 <- mat
mat2 <- mat2[apply(mat2, 2, function(y) !all(is.na(y))),]
mat2 <- mat2[,colSums(is.na(mat2))<nrow(mat2)]
dim(mat2)

ht2 = Heatmap(mat2,  col = f1, name = "",  row_title = "", column_title = paste0(dataset, ": GO sem sim score" ), show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , row_names_side ="left") # width = unit(140, "cm"),
ht2

pdf(paste( Project, "GOSemSim___ComplexHeatmap_to_T0_cluster", "Clust_purple.pdf", sep="_"), onefile=FALSE, width=12, height=10)
par(bg=NA)
draw(ht2, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

colnames(mat)[!colnames(mat) %in% colnames(mat2)]
# "1.vs.8"  "5.vs.8"  "1.vs.10" "8.vs.10" "0.vs.15" "6.vs.16" "8.vs.16"
dim(mat2)
dim(mat)
