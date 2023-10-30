#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()
options(bitmapType='cairo')


suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library("Seurat") 
  library(grDevices)
  #library(SCopeLoomR)
  options(bitmapType='cairo')
  library(dplyr)
  library(RColorBrewer)
})

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
Project <- "CTR_dscj1_harmony_R_I_sep"


message("--------------------------------------------------------------------------------")
message("+                    Load in sample table and matrix.su                         ")
message("+-------------------------------------------------------------------------------")


sampleTable <- read.csv("sampleTable_MERGED_53.csv")


Project <- "CTR_dscj1_harmony_R_500g"
matrix.su_sep <- readRDS("HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_R.Rds")
dataset_name <- "Remove_29_500g"
matrix.umap <- readRDS("Harmony_Remove_29_500g_matrix.umap_.Rds")

Project <- "CTR_dscj1_harmony_I_500g"
matrix.su_sep <- readRDS("HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_I.Rds")
dataset_name <- "Inhibit_32_500g"
matrix.umap <- readRDS("Harmony_Inhibit_32_500g_matrix.umap_.Rds")



matrix.su_T0 <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/matrix.su_T0_500g_harmony_3432.Rds")





Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.0.8
# I :::
#Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels =  c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0"))
# R :::
#Idents(matrix.su_sep) <- factor(Idents(matrix.su_sep), levels =   c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8"))

DotPlot(matrix.su_sep, features = c("Phf8","Elf2", "Elf3","Ascl2", "Gata1","Bag3", "Myc","Zfp644","Msx2","Clip1","Phactr4","Car2","Gjb2","E2f8","Gcm1","Mtmr7","Klf6","Cebpg","Npdc1","Ndrg1","Phlda2","Yy1","Tsix","Xist","Eed"), cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # ,"Ddx3y","Uty","Kdm5d" 


DotPlot(matrix.su_sep, features = c("Phf8","Elf2", "Elf3","Ascl2", "Gata1","Bag3", "Myc","Zfp644","Msx2","Clip1","Phactr4","Car2","Gjb2","E2f8","Gcm1","Mtmr7","Klf6","Cebpg","Npdc1","Ndrg1","Phlda2","Yy1","Tsix","Xist","Eed"), cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # ,"Ddx3y","Uty","Kdm5d" 



FeaturePlot(matrix.su_sep, features = c("Xist"))



# TEMPORA:::
cluster_order_I <- c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0")
cluster_order_R <- c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8")
# MONOCLE:::
#cluster_order_I <- c("9", "11" ,"16", "5" , "1" , "15", "6" , "12", "2" , "8" , "0" , "7" , "10", "14", "13", "4" , "3")
#cluster_order_R <- c("14", "13", "2" , "3" , "7", "0", "15", "6" , "12", "11", "5" , "10", "9" , "4" ,"8" ,"1")









dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")
Idents(dropseq.integrated) <- dropseq.integrated@meta.data$SCT_snn_res.0.8 
clusts_below100 <- as.data.frame(table(dropseq.integrated@meta.data[,"SCT_snn_res.0.8"]))
keep_clusters_over_100 <- clusts_below100[clusts_below100$Freq >= 100,]$Var1
dropseq.integrated <- subset(dropseq.integrated, SCT_snn_res.0.8 %in% keep_clusters_over_100)
dropseq.integrated@meta.data$SCT_snn_res.0.8 <- factor(dropseq.integrated@meta.data$SCT_snn_res.0.8 , levels = c("12", "15", "0", "3", "17", "11", "13", "2", "14","18", "7", "8", "5", "9", "10", "16", "1", "19", "6", "4")) # 14...19! # tempora





ALL_markers_T0 <- FindAllMarkers(matrix.su_T0)
write.csv(ALL_markers_T0, "CTR_dscj1_harmony_T0_500g_FindAllMarkers_wilcox_l2fc0.25_3432.csv")

#ALL_markers_I <- FindAllMarkers(matrix.su_sep)
#write.csv(ALL_markers_I, "CTR_dscj1_harmony_I_500g_FindAllMarkers_wilcox_l2fc0.25.csv")

ALL_markers_R <- FindAllMarkers(matrix.su_sep)
write.csv(ALL_markers_R, "CTR_dscj1_harmony_R_500g_FindAllMarkers_wilcox_l2fc0.25.csv")



message("--------------------------------------------------------------------------------")
message("+                CUSTOM Plot expression of selected genes on umap               ")
message("+-------------------------------------------------------------------------------")


PrepareSeuratForPlotting <- function(matrix.su, reso){
  head(matrix.su@meta.data)
  matrix.umap            <- as.data.frame(Embeddings(object=matrix.su, reduction="umap"))
  matrix.umap$Sample_2   <- rownames(matrix.umap) 
  matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
  sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
  rownames(sampleTable)  <- sampleTable$Sample_2 
  matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
  matrix.umap$Age        <- sampleTable[matrix.umap$Sample_2, ]$Age
  matrix.umap$Treatment  <- sampleTable[matrix.umap$Sample_2, ]$Treatment
  matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
  matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Batch
  matrix.umap$Project    <- sampleTable[matrix.umap$Sample_2, ]$Project
  
  clust               <- matrix.su@meta.data %>% dplyr::select((reso))  
  colnames(clust)     <- c("Cluster")
  matrix.umap$Cluster <- clust$Cluster
  head(matrix.umap)
  matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 
  unique(matrix.umap$Cluster)
  Integrated_data <- GetAssayData(matrix.su, assay.type = "SCT")
  return(list(matrix.umap, Integrated_data))
}




dropseq.integrated_ls <- PrepareSeuratForPlotting(dropseq.integrated, "SCT_snn_res.0.8")
head(dropseq.integrated_ls[[1]])

matrix.su_sep_ls <- PrepareSeuratForPlotting(matrix.su_sep, "SCT_snn_res.0.8")

matrix.su_T0_ls <- PrepareSeuratForPlotting(matrix.su_T0, "SCT_snn_res.0.4")



MakeCustomFeaturePlot <- function(gene2plot, matrix.su__ls, dataset_name, outdir){
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep="")) 
  gene_name <- print(gene2plot)
  matrix.umap <- matrix.su__ls[[1]]
  matrix.umap$selected_gene <- matrix.su__ls[[2]][gene2plot,]
  #matrix.umap$selected_gene <- GetAssayData(matrix.su, assay.type = "SCT")[gene2plot,]
  matrix.umap$selected_gene <- as.numeric(as.character(matrix.umap$selected_gene))
  matrix.umap <- matrix.umap[order(matrix.umap$selected_gene),]
  
  if(dataset_name == "Together"){
    matrix.umap$UMAP_1 <- matrix.umap$UMAP_1 * -1
  }
  
  myPalette <- colorRampPalette((brewer.pal(9, "YlGnBu")), bias = 1)
  sc <- scale_colour_gradientn(colours = myPalette(100) , limits=c(min(matrix.umap$selected_gene), max(matrix.umap$selected_gene))) 
  
  if(nrow(matrix.umap) > 20000){
  umap.gene    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point( aes(colour=selected_gene), alpha=0.3, size=0.3) +
    geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), linetype='dashed', colour="black", bins=3, alpha=0.5, size = 0.3) +
    sc +  coord_fixed() + xlab("") + ylab("") +  
    theme(title=element_text(size=6), text=element_text(size=6), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6)) +  theme_cowplot(12)  + theme(legend.title = element_blank())
  } else {
    umap.gene    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2)) +
      geom_point( aes(colour=selected_gene), alpha=0.7, size=0.9) +
      geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), linetype='dashed', colour="black", bins=3, alpha=0.5, size = 0.3) +
      sc +  coord_fixed() + xlab("") + ylab("") +  
      theme(title=element_text(size=6), text=element_text(size=6), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6)) +  theme_cowplot(12)  + theme(legend.title = element_blank())
  }

   return(umap.gene)
}


#umap_I_Phf8 <- MakeCustomFeaturePlot("Phf8", matrix.su_sep_ls, "Inhibit")
#umap_I_Gcm1 <- MakeCustomFeaturePlot("Gcm1", matrix.su_sep_ls, "Inhibit")
#umap_I_E2f8 <- MakeCustomFeaturePlot("E2f8", matrix.su_sep_ls, "Inhibit")
#umap_R_Phf8 <- MakeCustomFeaturePlot("Phf8", matrix.su_sep_ls, "Remove")
#umap_R_Gcm1 <- MakeCustomFeaturePlot("Gcm1", matrix.su_sep_ls, "Remove")
#umap_R_E2f8 <- MakeCustomFeaturePlot("E2f8", matrix.su_sep_ls, "Remove")

#pdf(paste( "CUSTOM_featurePlot_Harmony_500g_I_vs_R___LP_markers_Gcm1_E2f8_Phf8", ".pdf", sep=""), width=11, height=6)
#par(bg=NA)
#plot_grid(umap_I_Gcm1, umap_I_E2f8 , umap_I_Phf8, umap_R_Gcm1 , umap_R_E2f8, umap_R_Phf8, labels = c("Gcm1", "E2f8", "Phf8"), align = "hv", ncol = 3, label_size = 16, label_fontface = "plain")
#dev.off()



umap_orig_Elf2 <- MakeCustomFeaturePlot("Elf2", dropseq.integrated_ls, "Together")
Idents(dropseq.integrated) <- dropseq.integrated$Experiment
DotPlot(dropseq.integrated, features = c("Elf2", "Elf3", "Myc"), cols = c("green", "purple"))


# for plot in T0 fugure panel:::  Krt18, Rhox6, Rhox9 and Lgals3
umap_I_Krt18 <- MakeCustomFeaturePlot("Krt8", matrix.su_sep_ls, "Inhibit")
umap_T0_Krt18 <- MakeCustomFeaturePlot("Krt8", matrix.su_T0_ls, "T0")

umap_I_Rhox6 <- MakeCustomFeaturePlot("Rhox6", matrix.su_sep_ls, "Inhibit")
umap_T0_Rhox6 <- MakeCustomFeaturePlot("Rhox6", matrix.su_T0_ls, "T0")

umap_I_Rhox9 <- MakeCustomFeaturePlot("Rhox9", matrix.su_sep_ls, "Inhibit")
umap_T0_Rhox9 <- MakeCustomFeaturePlot("Rhox9", matrix.su_T0_ls, "T0")

umap_I_Lgals3 <- MakeCustomFeaturePlot("Lgals3", matrix.su_sep_ls, "Inhibit")
umap_T0_Lgals3 <- MakeCustomFeaturePlot("Lgals3", matrix.su_T0_ls, "T0")




plot_grid(umap_T0_Krt18, umap_T0_Rhox6, umap_T0_Rhox9, umap_T0_Lgals3 ,
          umap_I_Krt18, umap_I_Rhox6 , umap_I_Rhox9 , umap_I_Lgals3 , align = "hv", ncol = 4)

pdf(paste( "CUSTOM_featurePlot_Harmony_500g_I_vs_T0_GRID_2x4", "_genes", "_Krt18_Rhox6_Rhox9_Lgals3", "_3432.pdf", sep=""), width=12, height=6)
par(bg=NA)
plot_grid(umap_T0_Krt18, umap_T0_Rhox6, umap_T0_Lgals3 , umap_T0_Rhox9, 
          umap_I_Krt18, umap_I_Rhox6 ,  umap_I_Lgals3, umap_I_Rhox9 ,labels = c("Krt18", "Rhox6", "Rhox9","Lgals3"), align = "hv", ncol = 4, label_size = 16, label_fontface = "plain")
dev.off()










# FeaturePlot: Phf8, Gcm1, E2f8, ...
pdf(paste( "FeaturePlot", dataset, "_green_Phf8", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
FeaturePlot(matrix.su_sep,  features =c("Phf8","Gcm1","E2f8"), slot = "data", ncol = 3, pt.size =0.3 , max.cutoff  =  1.5)
dev.off()










