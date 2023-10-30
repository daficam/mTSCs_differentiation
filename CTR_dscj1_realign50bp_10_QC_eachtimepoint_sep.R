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

options(future.globals.maxSize = 4000 * 1024^2)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
Project <- "CTR_dscj1_MergedRealigned52_QC_sep_expt"

sampleTable <- read.csv("sampleTable_MERGED_53.csv")

message("--------------------------------------------------------------------------------")
message("+                                Load objects:::                                ")
message("+-------------------------------------------------------------------------------")

matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53.Rds")
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N721",])
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N729",])
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N729")
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N721")
head(matrix.su@meta.data,2)
matrix.su_R.1  <- subset(matrix.su, Experiment == "R.1")
matrix.su_R.4  <- subset(matrix.su, Experiment == "R.4")
matrix.su_R.24 <- subset(matrix.su, Experiment == "R.24")
matrix.su_R.36 <- subset(matrix.su, Experiment == "R.36")
matrix.su_R.48 <- subset(matrix.su, Experiment == "R.48")
matrix.su_I.1  <- subset(matrix.su, Experiment == "I.1")
matrix.su_I.4  <- subset(matrix.su, Experiment == "I.4")
matrix.su_I.24 <- subset(matrix.su, Experiment == "I.24")
matrix.su_I.36 <- subset(matrix.su, Experiment == "I.36")
matrix.su_I.48 <- subset(matrix.su, Experiment == "I.48")


matrix.su_list <- list(matrix.su_R.1 =matrix.su_R.1,
                       matrix.su_R.4 =matrix.su_R.4,
                       matrix.su_R.24=matrix.su_R.24,
                       matrix.su_R.36=matrix.su_R.36,
                       matrix.su_R.48=matrix.su_R.48,
                       matrix.su_I.1 =matrix.su_I.1,
                       matrix.su_I.4 =matrix.su_I.4,
                       matrix.su_I.24=matrix.su_I.24,
                       matrix.su_I.36=matrix.su_I.36,
                       matrix.su_I.48=matrix.su_I.48)


for(i in seq_along(matrix.su_list)){
  print(names(matrix.su_list)[i])
  matrix.su_list[[i]] <- PercentageFeatureSet(matrix.su_list[[i]], pattern = "^mt-", col.name = "percent.mt")
  matrix.su_list[[i]] <- subset(matrix.su_list[[i]], subset = percent.mt < 5)
  matrix.su_list[[i]] <- SCTransform(matrix.su_list[[i]], vars.to.regress = "percent.mt", verbose = TRUE)
  matrix.su_list[[i]] <- RunPCA(matrix.su_list[[i]], assay = "SCT", npcs = 30)
  matrix.su_list[[i]] <- RunUMAP(matrix.su_list[[i]], reduction = "pca", assay = "SCT", dims = 1:30, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 
  matrix.su_list[[i]] <- FindNeighbors(matrix.su_list[[i]], dims = 1:30)
  matrix.su_list[[i]] <- FindClusters(matrix.su_list[[i]], resolution = 0.2, n.start = 100, n.iter = 10, random.seed = 1) 
  gc()
}


#saveRDS(matrix.su_list[[1]], "matrix.su_R.1_QC.Rds")
#saveRDS(matrix.su_list[[2]], "matrix.su_R.4_QC.Rds")
#saveRDS(matrix.su_list[[3]], "matrix.su_R.24_QC.Rds")
#saveRDS(matrix.su_list[[4]], "matrix.su_R.36_QC.Rds")
#saveRDS(matrix.su_list[[5]], "matrix.su_R.48_QC.Rds")
#saveRDS(matrix.su_list[[6]], "matrix.su_I.1_QC.Rds")
#saveRDS(matrix.su_list[[7]], "matrix.su_I.4_QC.Rds")
#saveRDS(matrix.su_list[[8]], "matrix.su_I.24_QC.Rds")
#saveRDS(matrix.su_list[[9]], "matrix.su_I.36_QC.Rds")
#saveRDS(matrix.su_list[[10]], "matrix.su_I.48_QC.Rds")
#saveRDS(matrix.su_list, "matrix.su_list_qc_Epxt.Rds")

matrix.su_R.1  <- readRDS( "matrix.su_R.1_QC.Rds")
matrix.su_R.4  <- readRDS( "matrix.su_R.4_QC.Rds")
matrix.su_R.24 <- readRDS( "matrix.su_R.24_QC.Rds")
matrix.su_R.36 <- readRDS( "matrix.su_R.36_QC.Rds")
matrix.su_R.48 <- readRDS( "matrix.su_R.48_QC.Rds")
matrix.su_I.1  <- readRDS( "matrix.su_I.1_QC.Rds")
matrix.su_I.4  <- readRDS( "matrix.su_I.4_QC.Rds")
matrix.su_I.24 <- readRDS( "matrix.su_I.24_QC.Rds")
matrix.su_I.36 <- readRDS( "matrix.su_I.36_QC.Rds")
matrix.su_I.48 <- readRDS( "matrix.su_I.48_QC.Rds")
matrix.su_0.0  <- readRDS( "CTR_dscj1_realign52_matrix.su_T0_HET.Rds")



matrix.su_list <- readRDS("matrix.su_list_qc_Epxt.Rds")

sampleTable <- read.csv("sampleTable_MERGED_53.csv")

message("--------------------------------------------------------------------------------")
message("+                             Plot UMAP                                         ")
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


for(i in seq_along(matrix.su_list)){
matrix.su_list[[i]] <- CellCycleScoring(matrix.su_list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
}

names(matrix.su_list)
umap.batch_list <- list()
umap.Cluster_list <- list()
umap.cc_list <- list()
umap.batchSample_list <- list()

names(matrix.su_list)
#i <- 10
for(i in seq_along(matrix.su_list)){
  dataset_name <- names(matrix.su_list)[[i]]
  print(dataset_name)
  matrix.umap            <- as.data.frame(Embeddings(object=matrix.su_list[[i]], reduction="umap"))
  matrix.umap$Sample_2   <- rownames(matrix.umap) 
  matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","", matrix.umap$Sample_2)
  sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
  rownames(sampleTable)  <- sampleTable$Sample_2 
  matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
  matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
  matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
  matrix.umap$Project    <- sampleTable[matrix.umap$Sample_2, ]$Project
  #head(matrix.umap)

  cell_cycle               <- matrix.su_list[[i]]@meta.data %>% dplyr::select((Phase))  # 
  colnames(cell_cycle)     <- c("Phase")
  matrix.umap$Phase        <- cell_cycle$Phase

  reso                <- "res_0.2"
  clust               <- matrix.su_list[[i]]@meta.data %>% dplyr::select((SCT_snn_res.0.2))  # 
  colnames(clust)     <- c("Cluster")
  matrix.umap$Cluster <- clust$Cluster
  matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 
  #print(head(matrix.umap))
  unique(matrix.umap$Cluster)
  unique(matrix.umap$Batch)
  unique(matrix.umap$Phase)
  nrow(matrix.umap)
  unique(matrix.umap$Experiment)
  print(unique(matrix.umap$Batch))
  
  message("----------------------umap.Batch-sample-----------------------")
  
  matrix.umap$BatchSample <- paste(matrix.umap$Batch, matrix.umap$Sample, sep = "_")
  unique(matrix.umap$BatchSample)
  ncolours <- length(unique(matrix.umap$BatchSample))
  colours <- c("green", "red", "blue", "yellow","brown", "black", "magenta", "orange")
  colours <- colours[c(1:ncolours)]
  
  umap.BatchSample        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=BatchSample)) + geom_point(alpha=0.3, size=0.1) +
    ggtitle(paste0(dataset_name,  " BatchSample" )) +  coord_fixed() + xlab("") + ylab("") +
    scale_colour_manual("", values = colours) +
    theme(title=element_text(size=10), text=element_text(size=10),  axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
    guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + theme_classic()
  umap.batchSample_list[[i]] <- umap.BatchSample
  umap.BatchSample
  
message("----------------------umap.Cluster-----------------------")

umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.4, size=0.3) +
  #geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
  #                linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0(dataset_name, " QC :::  ", reso, " Cluster" )) + coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() 
umap.Cluster_list[[i]] <- umap.Cluster


message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.5) +
  #geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
  #                linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " Cell Cycle" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + theme_classic()
umap.cc_list[[i]] <- umap.cc


message("----------------------umap.Batch-----------------------")

#matrix.umap$Batch <- factor(matrix.umap$Batch, levels = c("AB0", "A00",  "ABC"))
ncolours <- length(unique(matrix.umap$Batch))
colours <- c("green", "red", "blue", "yellow","brown", "black", "magenta", "orange")
colours <- colours[c(1:ncolours)]

umap.Batch        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Batch)) + geom_point(alpha=0.3, size=0.1) +
  #geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
  #                linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0(dataset_name,  " Batch" )) +  coord_fixed() + xlab("") + ylab("") +
  #scale_colour_manual("", values = colours) +
  theme(title=element_text(size=10), text=element_text(size=10),  axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + theme_classic()
umap.batch_list[[i]] <- umap.Batch

}


#png(paste("ComplexHeatmap",  dataset_name,  "QC_batch.pdf", sep="_"),  width=1000, height=2000, type = "cairo") 
#par(bg=NA)
#plot_grid(umap.Batch, umap.cc, umap.Cluster, ncol = 1)
#dev.off()

plot_grid(umap.batch_list[[1]], umap.cc_list[[1]], umap.Cluster_list[[1]], umap.batchSample_list[[1]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[2]], umap.cc_list[[2]], umap.Cluster_list[[2]], umap.batchSample_list[[2]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[3]], umap.cc_list[[3]], umap.Cluster_list[[3]], umap.batchSample_list[[3]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[4]], umap.cc_list[[4]], umap.Cluster_list[[4]], umap.batchSample_list[[4]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[5]], umap.cc_list[[5]], umap.Cluster_list[[5]], umap.batchSample_list[[5]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[6]], umap.cc_list[[6]], umap.Cluster_list[[6]], umap.batchSample_list[[6]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[7]], umap.cc_list[[7]], umap.Cluster_list[[7]], umap.batchSample_list[[7]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[8]], umap.cc_list[[8]], umap.Cluster_list[[8]], umap.batchSample_list[[8]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[9]], umap.cc_list[[9]], umap.Cluster_list[[9]], umap.batchSample_list[[9]],  ncol = 2, vjust = 1, hjust = 1)
plot_grid(umap.batch_list[[10]], umap.cc_list[[10]], umap.Cluster_list[[10]], umap.batchSample_list[[10]],  ncol = 2, vjust = 1, hjust = 1)


plot_grid(umap.batch_list[[1]], umap.batch_list[[2]], umap.batch_list[[3]], umap.batch_list[[4]], umap.batch_list[[5]], umap.batch_list[[6]], umap.batch_list[[7]], umap.batch_list[[8]], umap.batch_list[[9]], umap.batch_list[[10]] , ncol = 3)

plot_grid(umap.batch_list[[1]], umap.batch_list[[2]], umap.batch_list[[3]], umap.batch_list[[6]], umap.batch_list[[7]], umap.batch_list[[8]] , ncol = 2, vjust = 1, hjust = 1)

plot_grid(umap.batch_list[[4]], umap.batch_list[[5]], umap.batch_list[[9]], umap.batch_list[[10]] , ncol = 2,vjust = 1, hjust = 1)

png(paste("ComplexHeatmap",  "sep_Expt",  "QC_batch.png", sep="_"),  width=1000, height=2000, type = "cairo") 
par(bg=NA)
print({plot_grid(umap.batch_list[[3]], umap.cc_list[[3]], umap.Cluster_list[[3]], ncol = 1)})
dev.off()



for(i in seq_along(matrix.su_list)){
  dataset_name <- names(matrix.su_list)[[i]]
  png(paste("ComplexHeatmap",  dataset_name,  "QC_batch.png", sep="_"),  width=1000, height=2000, type = "cairo") 
  par(bg=NA)
  print({plot_grid(umap.batch_list[[3]], umap.cc_list[[3]], umap.Cluster_list[[3]], ncol = 1)})
  dev.off()
}






message("--------------------------------------------------------------------------------")
message("+                    Detectable_genes...........                                ")
message("+-------------------------------------------------------------------------------")

matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53.Rds")
counts <- GetAssayData(matrix.su)
nrow(counts) # 20979

counts1 <- as.data.frame(counts[c(1:7000),])
counts2 <- as.data.frame(counts[c(7001:14000),])
counts3 <- as.data.frame(counts[c(14001:20979),])

filtered.counts1 <- counts1[rowSums(counts1>0)>=5, ]
dim(filtered.counts1) # 6966
dim(counts1)          # 7000 
filtered.counts2 <- counts2[rowSums(counts2>0)>=5, ]
dim(filtered.counts2) # 6980
dim(counts2)          # 7000 
filtered.counts3 <- counts3[rowSums(counts3>0)>=5, ]
dim(filtered.counts3) # 5841
dim(counts3)          # 7000 

expressed_genes <- c(rownames(filtered.counts1), rownames(filtered.counts2), rownames(filtered.counts3)) # 19787 genes
#saveRDS(expressed_genes, "background_expressed_genes_min_1count_5cells.Rds")

