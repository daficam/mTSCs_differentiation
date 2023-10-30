
rm(list=ls())
gc()
options(bitmapType='cairo')

######    module load R/4.0.2


library(BUSseq)
options(future.globals.maxSize = 4000 * 1024^2)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
resDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/BUSseq"
setwd(baseDir)
Project <- "CTR_dscj1_BUSseq"




message("--------------------------------------------------------------------------------")
message("+                  import data in right format::                                ")
message("+-------------------------------------------------------------------------------")

##### import data in right format::
#matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")
#matrix.su.list <- SplitObject(matrix.su, split.by = "Batch")
#expr_mat_all <- GetAssayData(matrix.su)
#dim(expr_mat_all) # 19755 217552
#expr_mat_list <- list(GetAssayData(matrix.su.list[[1]]), GetAssayData(matrix.su.list[[2]]),GetAssayData(matrix.su.list[[3]]),GetAssayData(matrix.su.list[[4]]),GetAssayData(matrix.su.list[[5]]),GetAssayData(matrix.su.list[[6]]),GetAssayData(matrix.su.list[[7]]))
#names(expr_mat_list) <- names(matrix.su.list)
#saveRDS(expr_mat_list, "BUSseq_expr_mat_list.Rds")

#var_genes <- VariableFeatures(matrix.su)
expr_mat_list_red <- list(as.matrix(GetAssayData(matrix.su.list[[1]][var_genes,])), 
                          as.matrix(GetAssayData(matrix.su.list[[2]][var_genes,])),
                          as.matrix(GetAssayData(matrix.su.list[[3]][var_genes,])),
                          as.matrix(GetAssayData(matrix.su.list[[4]][var_genes,])),
                          as.matrix(GetAssayData(matrix.su.list[[5]][var_genes,])),
                          as.matrix(GetAssayData(matrix.su.list[[6]][var_genes,])),
                          as.matrix(GetAssayData(matrix.su.list[[7]][var_genes,])))
names(expr_mat_list_red) <- names(matrix.su.list)
#saveRDS(expr_mat_list_red, "BUSseq_expr_mat_list_onlyVarGenes.Rds")



#Input data is a list
#CountData <- readRDS("BUSseq_expr_mat_list.Rds")
CountData <- readRDS("BUSseq_expr_mat_list_onlyVarGenes.Rds")
class(CountData) # [1] "list"
#The length of the list is three, so we have three batches
length(CountData) # 7
class(CountData[[1]])




message("--------------------------------------------------------------------------------")
message("+                            Model Fitting::                                    ")
message("+-------------------------------------------------------------------------------")

n.celltypes <- 22

#Conduct MCMC sampling and posterior inference for BUSseq model
BUSseqfits_res <- BUSseq_MCMC(ObservedData = CountData, seed = 1234, n.celltypes = n.celltypes , n.iterations = 500)

gc()
#saveRDS(BUSseqfits_res, paste0("BUSseqfits_res_n.celltypes_", n.celltypes,".Rds"))




message("--------------------------------------------------------------------------------")
message("+                     Extract results                                           ")
message("+-------------------------------------------------------------------------------")

setwd(resDir)
n.celltypes <- 20

BUSseqfits_res <- readRDS("BUSseqfits_res_n.celltypes_20.Rds")
#BUSseqfits_res <- readRDS("BUSseqfits_res_n.celltypes_22.Rds")
#BUSseqfits_res <- readRDS("BUSseqfits_res_n.celltypes_25.Rds")


summary(BUSseqfits_res)

celltypes_est <- celltypes(BUSseqfits_res)
saveRDS(celltypes_est, paste0("BUSseq_celltyes_est_n.celltypes_", n.celltypes,".Rds"))

location_batch_effects_est <- location_batch_effects(BUSseqfits_res)
saveRDS(location_batch_effects_est, paste0("BUSseq_location_batch_effects_est_n.celltypes_", n.celltypes,".Rds"))

overdispersion_est <- overdispersions(BUSseqfits_res)
saveRDS(overdispersion_est, paste0("BUSseq_overdispersion_est_n.celltypes_", n.celltypes,".Rds"))

cell_effects_est <- cell_effect_values(BUSseqfits_res)
saveRDS(cell_effects_est, paste0("BUSseq_cell_effects_est_n.celltypes_", n.celltypes,".Rds"))

celltype_mean_expression_est <- celltype_mean_expression(BUSseqfits_res)
## Each row represents a gene, and each column corresponds to a cell type.
head(celltype_mean_expression_est)
saveRDS(celltype_mean_expression_est, paste0("BUSseq_celltype_mean_expression_est_n.celltypes_", n.celltypes,".Rds"))

celltype_effects_est <- celltype_effects(BUSseqfits_res)
saveRDS(celltype_effects_est, paste0("BUSseq_celltype_effects_est_n.celltypes_", n.celltypes,".Rds"))




# Intrinsic Gene Identification

#obtain the intrinsic gene indicators
intrinsic_gene_indices <- intrinsic_genes_BUSseq(BUSseqfits_res)
#The estimated FDR, the first 500 genes are known intrinsic
#genes in the simulation setting.
false_discovery_ind <- !(intrinsic_gene_indices %in% 1:650)
fdr_est <- sum(false_discovery_ind)/length(intrinsic_gene_indices)
fdr_est

#  Corrected Read Count Data and Visualization
corrected_countdata <- corrected_read_counts(BUSseqfits_res)
## The output format is a "CountData" object with length equal to the batch number.
## Each element of the object is the corrected read count matrix.
## In each matrix, each row represents a gene and each column correspods to a cell.
class(corrected_countdata)
## [1] "CountData"
summary(corrected_countdata)

saveRDS(corrected_countdata, paste0("BUSseq_corrected_countdata_n.celltypes_", n.celltypes,".Rds"))







message("--------------------------------------------------------------------------------")
message("+                   Create batch corrected matrix.su                            ")
message("+-------------------------------------------------------------------------------")

setwd(resDir)

CountData <- readRDS("BUSseq_expr_mat_list_onlyVarGenes.Rds")
class(CountData)
CountData[[2]][1:5,1:5]

n.celltypes <- 20
BUSseqfits_res      <- readRDS( paste0("BUSseqfits_res_n.celltypes_", n.celltypes,".Rds"))
celltypes_est       <- readRDS( paste0("BUSseq_celltyes_est_n.celltypes_", n.celltypes,".Rds"))
corrected_countdata <- readRDS( paste0("BUSseq_corrected_countdata_n.celltypes_", n.celltypes,".Rds"))
corrected_countdata[[1]][1:5,1:5]
celltypes_est[[1]][1:5]


for(i in seq_along(corrected_countdata)){
   rownames(corrected_countdata[[i]]) <- rownames(CountData[[i]]) 
   colnames(corrected_countdata[[i]]) <- colnames(CountData[[i]])
   names(celltypes_est[[i]]) <- colnames(CountData[[i]])
}

CountData[[3]][1:5,1:5]
corrected_countdata[[3]][1:5,1:5]

celltypes_est_merged <- c(celltypes_est[[1]], celltypes_est[[2]], celltypes_est[[3]], celltypes_est[[4]], celltypes_est[[5]], celltypes_est[[6]], celltypes_est[[7]])
celltypes_est_merged[1:5]


matrix.su <- CreateSeuratObject(counts = Matrix(corrected_countdata[[1]], sparse = TRUE ) )
merge.su  <- CreateSeuratObject(counts = Matrix(corrected_countdata[[2]], sparse = TRUE ) )
matrix.su <- merge(matrix.su, merge.su )
merge.su  <- CreateSeuratObject(counts = Matrix(corrected_countdata[[3]], sparse = TRUE ) )
matrix.su <- merge(matrix.su, merge.su )
merge.su  <- CreateSeuratObject(counts = Matrix(corrected_countdata[[4]], sparse = TRUE ) )
matrix.su <- merge(matrix.su, merge.su )
merge.su  <- CreateSeuratObject(counts = Matrix(corrected_countdata[[5]], sparse = TRUE ) )
matrix.su <- merge(matrix.su, merge.su )
merge.su  <- CreateSeuratObject(counts = Matrix(corrected_countdata[[6]], sparse = TRUE ) )
matrix.su <- merge(matrix.su, merge.su )
merge.su  <- CreateSeuratObject(counts = Matrix(corrected_countdata[[7]], sparse = TRUE ) )
matrix.su <- merge(matrix.su, merge.su )
matrix.su_BUSseq <- matrix.su

head(matrix.su_BUSseq@meta.data,2)
#matrix.su_BUSseq@meta.data$clust_BUSseq <-  celltypes_est_merged
matrix.su_BUSseq@meta.data$clust_BUSseq <-  celltypes_est_merged[match( rownames(matrix.su_BUSseq@meta.data), names(celltypes_est_merged))]
matrix.su_BUSseq@meta.data$Cluster <- matrix.su_BUSseq@meta.data$clust_BUSseq


message("+            BIC QC                ")

# best number of cluster if BIC smallest!
BUSseqfits_res$BIC
# 919988394 for n.celltypes = 20
BIC_values <- c(919988394, ?, 917752558)
n.celltypes_values <- c(20,22,25)
plot(n.celltypes_values, BIC_values, xlab="cell type number", ylab="BIC", main="BIC plot", type="b")




message("--------------------------------------------------------------------------------")
message("+              matrix.su processing in Seurat and visualisation                 ")
message("+-------------------------------------------------------------------------------")

matrix.su_BUSseq <- ScaleData(matrix.su_BUSseq)
matrix.su_BUSseq <- RunPCA(matrix.su_BUSseq, npcs = 50, features = rownames(matrix.su_BUSseq))
matrix.su_BUSseq <- RunUMAP(matrix.su_BUSseq, reduction = "pca", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 
matrix.su_BUSseq <- FindNeighbors(matrix.su_BUSseq,  dims = 1:20) 
matrix.su_BUSseq <- FindClusters(matrix.su_BUSseq, resolution = 1.25)

head(matrix.su_BUSseq@meta.data,2)
saveRDS(matrix.su_BUSseq, "matrix.su_BUSseq_n.celltypes20.Rds")


message("+-------------------------------------------------------------------------------")
message("+                   Add missing metedata!!!!                                    ")
message("+-------------------------------------------------------------------------------")

colnames(matrix.su_BUSseq@meta.data)

phase_assignment <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/CellCycle_phase_assignment.Rds")

matrix.su_BUSseq@meta.data$Phase      <- phase_assignment[match(rownames(matrix.su_BUSseq@meta.data), rownames(phase_assignment)),]$Phase
matrix.su_BUSseq@meta.data$G2M.Score  <- phase_assignment[match(rownames(matrix.su_BUSseq@meta.data), rownames(phase_assignment)),]$G2M.Score
matrix.su_BUSseq@meta.data$S.Score    <- phase_assignment[match(rownames(matrix.su_BUSseq@meta.data), rownames(phase_assignment)),]$S.Score
matrix.su_BUSseq@meta.data$percent.mt <- phase_assignment[match(rownames(matrix.su_BUSseq@meta.data), rownames(phase_assignment)),]$percent.mt
head(matrix.su_BUSseq@meta.data,4)

sampleTable <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/sampleTable_MERGED_53.csv")
matrix.su_BUSseq@meta.data$Sample <- gsub( "_[ACTGN]{12}", "", rownames(matrix.su_BUSseq@meta.data))
matrix.su_BUSseq@meta.data$Batch  <- sampleTable[match(matrix.su_BUSseq@meta.data$Sample, sampleTable$sampleLabels2),]$Merged_batch
matrix.su_BUSseq@meta.data$Experiment  <- sampleTable[match(matrix.su_BUSseq@meta.data$Sample, sampleTable$sampleLabels2),]$CellType



message("+-------------------------------------------------------------------------------")
message("+                        Plotting UMAP - custom                                 ")
message("+-------------------------------------------------------------------------------")

dataset_name <- paste0("BUSseq_n.celltypes_", n.celltypes)
head(matrix.su_BUSseq@meta.data,2)

matrix.umap            <- as.data.frame(Embeddings(object=matrix.su_BUSseq, reduction="umap"))
#matrix.umap            <- as.data.frame(Embeddings(object=matrix.su_BUSseq, reduction="tsne"))
matrix.umap$Sample_2   <- rownames(matrix.umap) 
matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
rownames(sampleTable)  <- sampleTable$Sample_2 
matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
matrix.umap$Project    <- sampleTable[matrix.umap$Sample_2, ]$Project
head(matrix.umap)

cell_cycle             <- matrix.su_BUSseq@meta.data %>% dplyr::select((Phase))  # 
colnames(cell_cycle)   <- c("Phase")
matrix.umap$Phase      <- cell_cycle$Phase

reso                <- "clust_BUSseq_n.celltypes20"
clust               <- matrix.su_BUSseq@meta.data %>% dplyr::select((clust_BUSseq))  # 

colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 

head(matrix.umap,2)
unique(matrix.umap$Batch)
unique(matrix.umap$Cluster)


table(matrix.umap$Experiment)
table(matrix.umap$Cluster)

message("----------------------Cluster annotation -----------------------")

Clusters <- as.numeric(as.character(unique(matrix.umap$Cluster)))

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






message("----------------------umap.Cluster-----------------------")
unique(matrix.umap$Cluster)

matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18" ,  "19", "20" )) # , "21", "22", "23", "24", "25", "26"
# res 1.25:::
#matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18" ,  "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29" ))

umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.2, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0(dataset_name, " QC :::  ", reso, " Cluster" )) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme_classic() +
  annotate('text', label = paste0('5'), x =12.0166158632004, y = -0.354596237194373, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-4.62695646727445, y = -0.610337952625586, color = 'black', size=8) + annotate('text', label = paste0('15'), x =11.9851918176377, y = -0.332506934654547, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-6.13981342757109, y = 1.03517534350173, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.96622562849882, y = -0.620462516796424, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-5.98274374449613, y = -0.563924530994727, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-5.98562527144316, y = -0.599453786861731, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-6.00087595427397, y = -0.486044625294043, color = 'black', size=8) + annotate('text', label = paste0('4'), x =-5.98238373244169, y = -0.473897854339911, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-6.00903845274809, y = -0.446113327991797, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-5.98064876044157, y = -0.592272976887061, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-5.99711108649137, y = -0.536083082210852, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-5.98512316191557, y = -0.527166823398902, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-5.97901487791899, y = -0.551890591633155, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-5.98538113081816, y = -0.44350920106156, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-5.9747605367935, y = -0.525651136886908, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-6.05160475218656, y = -0.567672232639624, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-5.97538376295927, y = -0.454945305836035, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-6.06156897986296, y = -0.498675922405554, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-5.9943237348831, y = -0.501525501262976, color = 'black', size=8)
umap.Cluster





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
  theme(title=element_text(size=6), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) #
umap.Experiment




message("----------------------umap.Batch-----------------------")

# ABC 0B0 AB0 A0C 0BC A00 00C

umap.Batch    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Batch)) +
  geom_point(alpha=0.2, size=0.3) +
  ggtitle(paste0("Resolution ", reso, " :: Batch " )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Batch, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  #scale_colour_manual("", values = c("0.0"="black", "I.1"="yellow", "R.1"="yellow2", 
  #                                   "I.4"="orangered", "R.4"="red", "I.24"="lightblue2", "R.24"="lightblue4",
  #                                   "I.36"="blue", "R.36"="darkblue","I.48"="purple", "R.48"="orchid"), 
  #                    limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  theme_classic() + coord_fixed() +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=6), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) #
umap.Batch







