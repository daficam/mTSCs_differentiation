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
#Project <- "CTR_dscj1_MergedRealigned51_30L_rm004_9617_N729_R_I_sep"
Project <- "CTR_dscj1_MergedRealigned51_30L_rm004_9617_N729"
#Project <- "CTR_dscj1_MergedRealigned52_30L"

#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate__51__SCTClust_UMAP_30L_rm004_9617_N729_2.Rds")

message("--------------------------------------------------------------------------------")
message("+                                Load objects:::                                ")
message("+-------------------------------------------------------------------------------")

##matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformBatchScaledOut.Rds")
##matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53.Rds")
##matrix.su <- readRDS("INITIAL_NON_INTEGRATED/CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransform.Rds")
#matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")
matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate__51__SCTClust_UMAP_30L_rm004_9617_N729_2.Rds")
# Batch0004_9617_N721 and Batch0004_9617_N729 removed thus 51 samples now!!!

#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate__51__SCTClust_UMAP_30L_rm004_9617_N729.Rds")


nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N721",])
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N729",])

sampleTable <- read.csv("sampleTable_MERGED_53.csv")

matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N729")
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N721")


# sampleTable_183 <- read.csv("sampleTable2.csv")
# Batch0014_9339_N714, Batch0008_9339_N712, Batch0014_9339_N712, Batch0008_9616_N714, Batch0008_9616_N720, Batch0006_7635_N703, Batch0008_9339_N723 that were merged had align rate < 50%
#df <- merged_samples_molten[merged_samples_molten$value %in% c("Batch0014_9339_N714", "Batch0008_9339_N712", "Batch0014_9339_N712", "Batch0008_9616_N714", "Batch0008_9616_N720", "Batch0006_7635_N703", "Batch0008_9339_N723"),c(1,4)]
#df <- df[order(df$Name_of_rep_set),]
#unique(df$Name_of_rep_set)

View(matrix.su@meta.data)

nrow(matrix.su@meta.data)             # 229793 
min(matrix.su@meta.data$nFeature_RNA) # 300


matrix.su@meta.data$Sample <- matrix.su@meta.data$orig.ident
head(matrix.su@meta.data,2)
unique(matrix.su@meta.data$Sample) # 53 (some removed because of alignment rate too low!)
unique(matrix.su@meta.data$Batch)  # 7 batches as expected

DefaultAssay( matrix.su ) #<- "RNA"

table(matrix.su@meta.data$Experiment)
table(matrix.su@meta.data[c("SCT_snn_res.1.25","Experiment")])


nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0004_9617_N721",])
matrix.su <- subset(matrix.su, orig.ident2 != "Batch0004_9617_N721")


matrix.su@meta.data$Experiment <- factor(matrix.su@meta.data$Experiment, levels = c( "0.0", "I.1", "R.1",  "I.4", "R.4", "I.24", "R.24","I.36", "R.36","I.48", "R.48"  ))
DimPlot(matrix.su, group.by = "Experiment", split.by = "Batch", cols = c("black", "yellow", "yellow3",  "orangered", "red","lightblue2", "lightblue3","blue", "darkblue","purple", "orchid" ), ncol = 3, pt.size = 0.1)

Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.1.25


message("--------------------------------------------------------------------------------")
message("+                  tables batches per clust                                     ")
message("+-------------------------------------------------------------------------------")

colnames(dropseq.integrated@meta.data)
tmp_df <- dropseq.integrated@meta.data[,c("orig.ident" ,"Experiment", "Batch",  "SCT_snn_res.1.25")]

#tmp_df <- matrix.su@meta.data[,c("orig.ident" ,"Experiment", "Batch",  "SCT_snn_res.1.25")]
tmp_df$Experiment <- gsub( "\\.1", "\\.01" , tmp_df$Experiment)
tmp_df$Experiment <- gsub( "\\.4", "\\.04" , tmp_df$Experiment)
tmp_df$Experiment <- gsub( "\\.048", "\\.48" , tmp_df$Experiment)
tmp_df$Expt_batch <- paste(tmp_df$Experiment, tmp_df$Batch, sep = "_")
head(tmp_df)
tmp_df$integrated_snn_res.1.25 <- as.numeric(as.character(tmp_df$SCT_snn_res.1.25))

tmp_df <- as.data.frame(table(tmp_df[,c("SCT_snn_res.1.25", "Batch")]))
tmp_df2 <- subset(tmp_df, tmp_df$Freq > 30)
tmp_df2 <- as.data.frame(table(tmp_df2$SCT_snn_res.1.25))
mean(tmp_df2$Freq)


message("--------------------------------------------------------------------------------")
message("+                  tables cells per clust / expt-batch                          ")
message("+-------------------------------------------------------------------------------")

tmp_df <- matrix.su@meta.data[,c("orig.ident" ,"Experiment", "Batch", "SCT_snn_res.0.8", "SCT_snn_res.1", "SCT_snn_res.1.25", "Phase")]
tmp_df$Experiment <- gsub( "\\.1", "\\.01" , tmp_df$Experiment)
tmp_df$Experiment <- gsub( "\\.4", "\\.04" , tmp_df$Experiment)
tmp_df$Experiment <- gsub( "\\.048", "\\.48" , tmp_df$Experiment)
tmp_df$Expt_batch <- paste(tmp_df$Experiment, tmp_df$Batch, sep = "_")
head(tmp_df)
tmp_df$SCT_snn_res.0.8 <- as.numeric(as.character(tmp_df$SCT_snn_res.0.8))
tmp_df$SCT_snn_res.1 <- as.numeric(as.character(tmp_df$SCT_snn_res.1))
tmp_df$SCT_snn_res.1.25 <- as.numeric(as.character(tmp_df$SCT_snn_res.1.25))
table(tmp_df[,c("Expt_batch","SCT_snn_res.0.8")])
table(tmp_df[,c("Expt_batch","SCT_snn_res.1")])
table(tmp_df[,c("Expt_batch","SCT_snn_res.1.25")])

message("--------------------------------------------------------------------------------")
message("+                  tables cells per clust / expt-batch                          ")
message("+-------------------------------------------------------------------------------")

tmp_df_t0 <- subset(tmp_df, tmp_df$Experiment == "0.0")
head(tmp_df_t0)
table(tmp_df_t0[,c("Phase","SCT_snn_res.1")])

tmp_phase <- as.data.frame(table(tmp_df_t0[,c("orig.ident","Phase")]))

ggplot(data=tmp_phase, aes(x=orig.ident, y=Freq, fill=Phase)) + geom_bar(stat="identity") +   theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


tmp_clust <- as.data.frame(table(tmp_df_t0[,c("orig.ident","SCT_snn_res.1")]))
tmp_clust <-tmp_clust[tmp_clust$Freq > 10, ]
ggplot(data=tmp_clust, aes(x=orig.ident, y=Freq, fill=SCT_snn_res.1)) + geom_bar(stat="identity") +   theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tmp_clust_phase <-  as.data.frame(table(tmp_df_t0[,c("Phase","SCT_snn_res.1")]))
ggplot(data=tmp_clust_phase, aes(x=SCT_snn_res.1, y=Freq, fill=Phase)) + geom_bar(stat="identity") +   theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




message("--------------------------------------------------------------------------------")
message("+                           sc to bulk                                          ")
message("+-------------------------------------------------------------------------------")

library(DESeq2)

#SPLIT <- "Batch"
SPLIT <- "Sample"


sc_to_bulk_obj <- SplitObject(matrix.su, split.by = SPLIT)
names(sc_to_bulk_obj)


list_bulk <- list()
for (i in seq_along(sc_to_bulk_obj)){
  list_bulk[[i]] <- rowSums(as.matrix(GetAssayData(sc_to_bulk_obj[[i]], slot = "counts")))
}
names(list_bulk) <- names(sc_to_bulk_obj)

df <- as.data.frame(t(data.frame(matrix(unlist(list_bulk), nrow=length(list_bulk), byrow=T))))
colnames(df) <- names(list_bulk)
rownames(df) <- names(list_bulk[[1]])
head(df,2)


sampleTable_sc2bulk <- sampleTable


if(SPLIT == "Batch"){
  sampleTable_sc2bulk <- sampleTable_sc2bulk[!duplicated(sampleTable_sc2bulk$Batch),]
  rownames(sampleTable_sc2bulk) <- sampleTable_sc2bulk$Batch
  rownames(sampleTable_sc2bulk) <- gsub( "Batch", "Dropseq_", rownames(sampleTable_sc2bulk))
} else if(SPLIT == "Sample"){
  sampleTable_sc2bulk <- sampleTable[!duplicated(sampleTable$sampleLabels),]
  rownames(sampleTable_sc2bulk) <- sampleTable_sc2bulk$sampleLabels
}


sampleTable_sc2bulk <- sampleTable_sc2bulk[match(colnames(df), rownames(sampleTable_sc2bulk)),]


dds <- DESeq2:: DESeqDataSetFromMatrix(countData=df, colData=sampleTable_sc2bulk, design=~CellType)
colData(dds)


#rld <- DESeq2::rlogTransformation(dds, blind=TRUE)
vst <- DESeq2::vst(dds, blind=TRUE)
#saveRDS(vst, "CTR_dscj1_Dropseq_newAlignment50bp_MERGED_53_sc-to-bulk__vst.Rds")



message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

#vst <- readRDS("CTR_dscj1_Dropseq_newAlignment50bp_MERGED_53_sc-to-bulk__vst.Rds")
rld <- vst

elementTextSize <- 14
pca = prcomp(t(assay(rld)))
rv = rowVars(assay(rld))


pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

sampleTable_sc2bulk$ID <- paste(sampleTable_sc2bulk$sampleLabels, sampleTable_sc2bulk$CellType, sep = "__")
scores <- data.frame(pca$x, sampleTable_sc2bulk)
#scores2 <- scores
#scores <- subset(scores2, scores2$uniq.align.rates > 30 )

Supp_Fig_1_A <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(Merged_batch) ))) + 
  geom_point(size = 5 , alpha = 0.5) + 
  xlab(pc1lab) + ylab(pc2lab) + 
  scale_colour_manual(name="Batch", values = c("green", "blue", "violet", "red", "orange", "yellow", "black"))  +
  theme_classic() 
  #geom_text_repel(aes(label=ID), col = "black") 
  #geom_text_repel(aes(label=ID), col = "black") 
Supp_Fig_1_A


scores$CellType <- factor(scores$CellType, levels = c("0.0","R.1", "R.4","R.24" , "R.36","R.48" , "I.1","I.4" , "I.24" ,"I.36", "I.48", "15980_xN727", "15987_xN715"  ))
#scores$CellType <- factor(scores$CellType, levels = c("0.0","R.1", "R.4","R.24" , "R.36","R.48" , "I.1","I.4" , "I.24" ,"I.36", "I.48"  ))

Supp_Fig_1_b <- ggplot(scores, aes(x = PC1, y = PC2, col = CellType )) + 
  geom_point(size = 5 , alpha = 0.7) + 
  xlab(pc1lab) + ylab(pc2lab) + 
  scale_colour_manual(name="Batch", values = c("0.0"="black", "I.1"="yellow", "R.1"="yellow3",  "I.4"="orangered", "R.4"="red", "I.24"="lightblue2", "R.24"="lightblue3","I.36"="blue", "R.36"="darkblue","I.48"="purple", "R.48"="orchid", "15980_xN727" = "grey85", "15987_xN715" = "grey" ))  +
  theme_classic() 
Supp_Fig_1_b




message("+-------------------------------------------------------------------------------")
message("+                         Calculate Correlation Matrix                          ")
message("+-------------------------------------------------------------------------------")

dds <- estimateSizeFactors(dds)
s2c <- sampleTable_sc2bulk

elementTextSize <- 10


TOPNUM <- 1000
rv     <- rowVars(as.matrix(assay(rld)))
select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
pca    <- prcomp(as.matrix(assay(rld))[select,])

norm_mat_sel <- as.matrix(assay(rld))[select,]
#norm_mat_sel <- cbind(group = col_ann[,2] , as.data.frame(t(scaled_matrix[select, ])))
#norm_mat_sel <- as.data.frame(t(scaled_matrix[select, ]))
#pca    <- prcomp(norm_mat_sel[,-1])

#  http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html


sampleTable <- sampleTable[order(sampleTable$CellType),]
sampleTable <- sampleTable[order(sampleTable$Age),]

norm_mat_sel2 <- norm_mat_sel[ , match( sampleTable$sampleLabels, colnames(norm_mat_sel) )]
colnames(norm_mat_sel2) <- paste(colnames(norm_mat_sel2), sampleTable$CellType, sep = "_")
norm_mat_sel2[1:5,1:5]

plot(hclust(as.dist(1-cor(norm_mat_sel))))

corr_mat <- cor(norm_mat_sel2)
min(corr_mat)

library(ggfortify)
ggplot2::autoplot(corr_mat) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



hmisc_mat <- Hmisc::rcorr(norm_mat_sel2, type=c("pearson"))		#Matrix of Correlations and P-values
hmisc_mat$P
hmisc_mat$r

library(circlize)
library(ComplexHeatmap)

f1 = colorRamp2( c(0, 0.7, 0.85, 1), c("white", "lightskyblue",  "deepskyblue3", "blue4"), space = "RGB") 

ht1 = Heatmap(as.matrix(hmisc_mat$r),  col = f1, name = "sc-to-Bulk",  row_title = "", column_title = "sc-to-Bulk corr", show_row_names = TRUE, heatmap_legend_param = list(title = "Correlation", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = T, cluster_rows = T , row_names_side ="left")
ht1

pdf(paste("ComplexHeatmap",  "MERGED_newAlignment50bp_CLUST", "Corr_sc_to_bulk.pdf", sep="_"), onefile=FALSE, width=30, height=30) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()






message("+-------------------------------------------------------------------------------")
message("+                          initial Seurat processing                            ")
message("+-------------------------------------------------------------------------------")

# store mitochondrial percentage in object meta data
matrix.su <- PercentageFeatureSet(matrix.su, pattern = "^mt-", col.name = "percent.mt")
max(matrix.su@meta.data$percent.mt) # 73.59477%
#View(matrix.su@meta.data)
dim(matrix.su) #  20979 229793


# Visualize QC metrics as a violin plot
Idents(matrix.su) <- matrix.su@meta.data$Batch
VlnPlot(matrix.su, features = c("nFeature_RNA", "percent.mt"), ncol = 2, pt.size = 0.05, assay = "RNA")
mean(matrix.su@meta.data$nFeature_RNA)
median(matrix.su@meta.data$nFeature_RNA)

# subset MT 
matrix.su <- subset(matrix.su, subset = percent.mt < 5)
dim(matrix.su) #  20979 217655


# run sctransform

#matrix.su <- NormalizeData(matrix.su)
#matrix.su <- FindVariableFeatures(matrix.su, nfeatures = 3000)
#matrix.su <- ScaleData(matrix.su, vars.to.regress = "percent.mt")

matrix.su <- SCTransform(matrix.su, return.only.var.genes = FALSE, vars.to.regress = c("percent.mt"), verbose = T) # , "Batch"
#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformBatchScaledOut.Rds")

# These are now standard steps in the Seurat workflow for visualization and clustering
matrix.su <- RunPCA(matrix.su, assay = "SCT", npcs = 50)

#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 25L) 
#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 28L) 
matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L) 
#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT" , dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.005, n.neighbors= 35L) 
#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT" , dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.005, n.neighbors= 40L) 
#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT" , dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.005, n.neighbors= 10L) 

#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")

#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_40L.Rds")

#matrix.su <- RunUMAP(matrix.su, dims = 1:50, min.dist = 0.01 )
#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_std_dim50.Rds")


#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L,  local.connectivity = 20L) # spread = 2,
#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L_colConn20L.Rds")

#matrix.su <- RunUMAP(matrix.su, reduction = "pca", assay = "SCT", dims = 1:50, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors= 30L,  spread = 2) # spread = 2,
#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate__51__SCTClust_UMAP_30L_rm004_9617_N729.Rds")


matrix.su <- RunTSNE(matrix.su, perplexity = 50, do.fast = TRUE)




#matrix.su <- FindClusters(matrix.su)



# some visualisation::

DimPlot(matrix.su, label = TRUE, reduction = "pca") + NoLegend()
VizDimLoadings(matrix.su, dims = 1:2, reduction = "pca")
print(matrix.su[["pca"]], dims = 1:5, nfeatures = 10)

Idents(matrix.su) <- matrix.su@meta.data$Experiment
#Idents(matrix.su) <- matrix.su@meta.data$Batch
DimPlot(matrix.su, label = TRUE) + NoLegend()

ElbowPlot(matrix.su)



message("--------------------------------------------------------------------------------")
message("+               cell cycle with SEURAT for T.0                                  ")
message("+-------------------------------------------------------------------------------")


#exp.mat <- read.table(file = "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/SEURAT_172_samples_0006_0007_0008_0011_0014/INTEGRATED_SEURAT_172/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)
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





message("+-------------------------------------------------------------------------------")
message("+                      Find Clusters and Clustree                               ")
message("+-------------------------------------------------------------------------------")

matrix.su <- FindNeighbors(matrix.su, dims = 1:50, force.recalc = TRUE)
matrix.su <- FindClusters(matrix.su, resolution = 0,    dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1)  
matrix.su <- FindClusters(matrix.su, resolution = 0.2,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
matrix.su <- FindClusters(matrix.su, resolution = 0.6,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
matrix.su <- FindClusters(matrix.su, resolution = 1,    dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
matrix.su <- FindClusters(matrix.su, resolution = 1.5,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
matrix.su <- FindClusters(matrix.su, resolution = 2,    dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
matrix.su <- FindClusters(matrix.su, resolution = 1.25,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 
matrix.su <- FindClusters(matrix.su, resolution = 0.8,  dims = 1:20 , n.start = 100, n.iter = 10, random.seed = 1) 

colnames(matrix.su@meta.data)



tree_to_plot <- clustree::clustree(matrix.su, prefix = "SCT_snn_res.") # "RNA_snn_res."

pdf(paste("matrix.su_MERGED_REALIGNED_ini_clustree.pdf", sep=""), width=12,height=8)
par(bg=NA)
tree_to_plot
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                        Plotting UMAP - custom                                 ")
message("+-------------------------------------------------------------------------------")

dataset_name <- "matrix.su_MERGED_REALIGNED_51_ini"
#dataset_name <- "matrix.su_MERGED_REALIGNED_52_ini"
head(matrix.su@meta.data,2)

matrix.umap            <- as.data.frame(Embeddings(object=matrix.su, reduction="umap"))
#matrix.umap            <- as.data.frame(Embeddings(object=matrix.su, reduction="tsne"))
matrix.umap$Sample_2   <- rownames(matrix.umap) 
matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
rownames(sampleTable)  <- sampleTable$Sample_2 
matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
matrix.umap$Project    <- sampleTable[matrix.umap$Sample_2, ]$Project
head(matrix.umap)

cell_cycle             <- matrix.su@meta.data %>% dplyr::select((Phase))  # 
colnames(cell_cycle)   <- c("Phase")
matrix.umap$Phase      <- cell_cycle$Phase

reso                <- "SCT_snn_res.1.25"
clust               <- matrix.su@meta.data %>% dplyr::select((SCT_snn_res.1.25))  # 

colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 

head(matrix.umap,2)
unique(matrix.umap$Batch)
unique(matrix.umap$Cluster)


table(matrix.umap$Experiment)
table(matrix.umap$Cluster)

message("----------------------QC small dots on UMAP-----------------------")

check_df <- subset(matrix.umap, matrix.umap$UMAP_1 < -5 & matrix.umap$UMAP_2 >2 )
check_df1 <- subset(check_df, check_df$UMAP_2 > 3)
check_df2 <- subset(check_df, check_df$UMAP_2 < 3)
table(check_df2[,c("Experiment", "Cluster")])
table(check_df1[,c("Experiment", "Cluster")])
table(check_df2$Cluster)
table(check_df2$Experiment)
table(check_df1$Cluster)
table(check_df1$Experiment)



message("----------------------Cluster annotation -----------------------")

Clusters <- as.numeric(levels(matrix.umap$Cluster))

coord_x <- list()
coord_y <- list()

for (i in 1:length(Clusters)){
  coord_x[i] <- median(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$UMAP_1)  
  coord_y[i] <- median(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$UMAP_2)  
  #coord_x[i] <- median(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$tSNE_1 )  
  #coord_y[i] <- median(matrix.umap[matrix.umap$Cluster ==  Clusters[i],]$tSNE_2 ) 
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
# res 0.8:::
# matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18" )) 
# res 1.0:::
matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18" ,  "19", "20", "21", "22", "23", "24", "25", "26","27","28","29","30","31" ))
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
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)


#umap.Cluster

#umap.Cluster_1 <- umap.Cluster
#umap.Cluster_0.8 <- umap.Cluster
umap.Cluster_1.25 <- umap.Cluster

png(paste0("UMAP30_MergedRealigned_52_ini_", Project, "_ggplot_", "Clusters_grid", "__30L_",".png"), width=2000, height=1200, type="cairo")
par(bg=NA)
plot_grid(umap.Cluster_0.8, umap.Cluster_1, umap.Cluster_1.25,umap.Experiment, ncol = 2, nrow = 2, align = "hv")
dev.off()


png(paste0("tSNE_MergedRealigned_52_ini_", Project, "_ggplot_", "plot_grid_clust_expt_cc",".png"), width=3000, height=1000, type="cairo")
par(bg=NA)
plot_grid(umap.Cluster, umap.Experiment, umap.cc,umap.Experiment, ncol = 3, nrow = 1, align = "hv")
dev.off()


png(paste0("UMAP30_MergedRealigned_52_ini_", Project, "_ggplot_", "Cluster", reso, "__30L_",".png"), width=1000, height=900, type="cairo")
#pdf(paste0("UMAP30_MergedRealigned_53_ini_", Project, "_ggplot_", "Expt", ".pdf"), width=10,height=9)
par(bg=NA)
umap.Cluster
dev.off()



message("----------------------umap.QC_1_sample-----------------------")

table(matrix.su@meta.data[,c(4,5)])

matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Dups_setA1_reps_N724", "Dups_setA1_reps_N724", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Dups_setA2_reps_N719", "Dups_setA2_reps_N719", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Dups_setB_reps_N704", "Dups_setB_reps_N704", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0004_9617_N718", "Batch0004_9617_N718", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0004_9617_N719", "Batch0004_9617_N719", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0004_9617_N720", "Batch0004_9617_N720", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0004_9617_N722", "Batch0004_9617_N722", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0004_9617_N726", "Batch0004_9617_N726", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0004_9617_N727", "Batch0004_9617_N727", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0004_9617_N729", "Batch0004_9617_N729", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Sample_2 == "Batch0011_15988_N704", "Batch0011_15988_N704", matrix.umap$QC)

#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc_0        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.5, size=0.3) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 samples/batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "Dups_setA1_reps_N724" = "blue1", "Dups_setA2_reps_N719"="lightblue2", "Dups_setB_reps_N704" = "lightblue4", "Batch0004_9617_N729" = "black", "Batch0004_9617_N727" = "green1", "Batch0004_9617_N726" = "brown", "Batch0004_9617_N722" = "limegreen", "Batch0004_9617_N720" = "violet","Batch0004_9617_N719" = "orange","Batch0004_9617_N718" = "red", "Batch0011_15988_N704" = "yellow")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic()  + xlim(c(6,7)) + ylim(c(-3, -1.8)) # xlim(c(3,9)) + ylim(c(-3.5, 2.5)) # 
umap.qc_0


message("----------------------umap.QC_1_sample-----------------------")
table(matrix.su@meta.data[,c( "SCT_snn_res.1.25", "Experiment")])

table(matrix.su@meta.data[,c(5,8)])

matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "A00", "A00_t0", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "00C", "00C_t0", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0B0", "0B0_t0", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0BC", "0BC_t0", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc_0        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.3, size=0.1) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "A00_t0" = "blue1", "00C_t0"="red", "0B0_t0" = "orange", "0BC_t0" = "green")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic()  + xlim(c(3,9)) + ylim(c(-3.5, 2.5)) # xlim(c(6,7)) + ylim(c(-3, -1.8)) #
umap.qc_0


matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1", "I.1", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "0B0", "0B0_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "A00", "A00_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "A0C", "A0C_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "ABC", "ABC_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "AB0", "AB0_I.1", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc_I.1        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC I.1 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "A00_I.1" = "blue1", "A0C_I.1"="red", "0B0_I.1" = "orange", "ABC_I.1" = "green", "AB0_I.1"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc_I.1




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24", "I.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$Batch == "ABC", "ABC_I.24", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$Batch == "AB0", "AB0_I.24", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc_I.24        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_I.24" = "green", "AB0_I.24"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc_I.24



matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4", "R.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "ABC", "ABC_R.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "AB0", "AB0_R.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "A0C", "A0C_R.4", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc_R.4        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_R.4" = "green", "AB0_R.4"="purple", "A0C_R.4"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc_R.4




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4", "I.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "ABC", "ABC_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "AB0", "AB0_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "A0C", "A0C_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "A00", "A00_I.4", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc_I.4        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_I.4" = "green", "AB0_I.4"="purple", "A0C_I.4"="blue", "A00_I.4" = "orange")) +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc_I.4




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24", "R.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24" & matrix.umap$Batch == "ABC", "ABC_R.24", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24" & matrix.umap$Batch == "A00", "A00_R.24", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc_R.24        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_R.24" = "green", "A00_R.24"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc_R.24




#png(paste0("UMAP30_MergedRealigned_53_ini_", Project, "_ggplot_", "CellCycle", ".png"), width=1000,height=900, type="cairo")
#par(bg=NA)
#umap.cc
#dev.off()




message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.1) +
  #geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                 # linetype='dashed', colour="black", bins=3, alpha=0.5) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " Cell Cycle" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() #   theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +

umap.cc

png(paste0("UMAP30_MergedRealigned_51_ini_", Project, "_ggplot_", "CellCycle","__30L_", ".png"), width=1000,height=900, type="cairo")
par(bg=NA)
umap.cc
dev.off()


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
  theme(title=element_text(size=6), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)

umap.Experiment

png(paste0("UMAP30_MergedRealigned_51_ini_", Project, "_ggplot_", "Expt", "__30L_",".png"), width=1000, height=900, type="cairo")
#pdf(paste0("UMAP30_MergedRealigned_53_ini_", Project, "_ggplot_", "Expt", ".pdf"), width=10,height=9)
par(bg=NA)
umap.Experiment
dev.off()








message("----------------------umap.t.0-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_timepoints")
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "A00", "A00_t0", matrix.umap$Mixer)
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "00C", "00C_t0", matrix.umap$Mixer)
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0B0", "0B0_t0", matrix.umap$Mixer)
matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$Batch == "0BC", "0BC_t0", matrix.umap$Mixer)
table(matrix.umap$Mixer)

umap.t0    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.2, size=0.2) +
  ggtitle(paste0("Resolution ", reso, " T0 in Batches" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "A00_t0" = "blue1", "00C_t0"="red", "0B0_t0" = "orange", "0BC_t0" = "green")) +
  coord_fixed() +
  xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)




message("----------------------umap.t.1-----------------------")

umap.Experiment1    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="blue", "R.1"="red", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)


message("----------------------umap.t.4-----------------------")

umap.Experiment4    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="blue", "R.4"="red", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)




message("----------------------umap.t.24-----------------------")

umap.Experiment24    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="blue", "R.24"="red",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c('0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))  + 
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)



message("----------------------umap.t.36-----------------------")

umap.Experiment36    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="blue", "R.36"="red","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)





message("----------------------umap.t.48-----------------------")

umap.Experiment48    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="blue", "R.48"="red"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =-3.50811649617225, y = 0.834185918023212, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.47401810940772, y = 0.85874511783324, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-0.732355726328196, y = 2.45153876130782, color = 'black', size=8) + annotate('text', label = paste0('11'), x =5.37852762881249, y = -0.894645522187131, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-0.716056478586496, y = 2.95860524957381, color = 'black', size=8) + annotate('text', label = paste0('13'), x =-1.38571680840522, y = 4.26471718137465, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.23336864766151, y = -3.62435392553605, color = 'black', size=8) + annotate('text', label = paste0('15'), x =0.450673418674169, y = -1.16296065861978, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-2.46087325867683, y = 2.99651261155806, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-0.475481999483408, y = -2.73982087785997, color = 'black', size=8) + annotate('text', label = paste0('18'), x =4.87011860552758, y = 0.394715626885517, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-3.28308297452003, y = -1.2503695588807, color = 'black', size=8) + annotate('text', label = paste0('2'), x =5.42463778201073, y = -1.93996957952775, color = 'black', size=8) + annotate('text', label = paste0('20'), x =-1.71511782464057, y = -1.06195120031633, color = 'black', size=8) + annotate('text', label = paste0('21'), x =6.35740231219262, y = 0.84672768895827, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-6.59089280423194, y = -0.693637425849812, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-5.06122089680702, y = -0.806639204571622, color = 'black', size=8) + annotate('text', label = paste0('24'), x =-0.855464947786631, y = -1.47162834818162, color = 'black', size=8) + annotate('text', label = paste0('25'), x =8.1973671788308, y = -3.04469076807298, color = 'black', size=8) + annotate('text', label = paste0('26'), x =0.187999385032354, y = 2.12382252996169, color = 'black', size=8) + annotate('text', label = paste0('27'), x =-6.80939866360694, y = 0.865519245316608, color = 'black', size=8) + annotate('text', label = paste0('28'), x =-0.144381879274251, y = 0.804538925339801, color = 'black', size=8) + annotate('text', label = paste0('29'), x =0.431094514760672, y = -2.33914385254182, color = 'black', size=8) + annotate('text', label = paste0('3'), x =0.21795560661286, y = 3.6548932391902, color = 'black', size=8) + annotate('text', label = paste0('30'), x =-6.95299388226539, y = 4.08046181981765, color = 'black', size=8) + annotate('text', label = paste0('31'), x =1.30061970892876, y = -2.7303709190587, color = 'black', size=8) + annotate('text', label = paste0('4'), x =2.35387180987328, y = -1.37333313638963, color = 'black', size=8) + annotate('text', label = paste0('5'), x =1.99907087031334, y = -2.27119843179979, color = 'black', size=8) + annotate('text', label = paste0('6'), x =0.692084717187582, y = -3.73428861315049, color = 'black', size=8) + annotate('text', label = paste0('7'), x =-1.76061869915992, y = 0.784748752762897, color = 'black', size=8) + annotate('text', label = paste0('8'), x =6.48919581118554, y = -0.990501384089367, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-4.13335276898414, y = -0.214664126584904, color = 'black', size=8)




png(paste( "UMAP75", Project, "ggplot", "Expt_grid_eachTimePoint", reso,"min.dist0.01_nneigh30L", ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
plot_grid(umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 2, align = "hv")
dev.off()







message("+-------------------------------------------------------------------------------")
message("+                 plot previously identified genes                              ")
message("+-------------------------------------------------------------------------------")


# top QC check genes: 
FeaturePlot(matrix.su,  reduction = "tsne", features = c("Elf5", "Eomes", "Sox2","Ascl2","Gcm1","Ovol2", "Anxa1" , "Plac1", "Uba6")) # "Anxa1,"Pparg" "Cdx2"




FeaturePlot(matrix.su, features = c("Elf5", "Cdx2","Eomes","Tead4","Sox2","Ascl2","Gcm1","Ovol2"))
# LP ("Gcm1","Ovol2"), EPC = Ascl2
FeaturePlot(matrix.su, features = c("Usf2", "Nfya", "Yy1", "Maff", "Gcm1","Ovol2", "Pparg","Ghrl1")) # LP
FeaturePlot(matrix.su, features = c("Vezf1", "Gata1", "Klf7", "Bdp1","Ascl2")) # EPC
FeaturePlot(matrix.su, features = c("Bhlhe40", "Smarca4", "Myc", "Foxo3","Sp3","Srebf2","Hdac6","Tead3", "Ppard")) # TFs
# other TFs : "E2f8","","Grhl1","","Foxo4","Hdac2","","Junb","Nfkb2","Creb3l2","Pou2f1","Mycn","","Max","","","E2f7","Fos","Elf2"


VlnPlot(matrix.su, features = c("Lgals3", "Krt18"), assay = "SCT", log = TRUE)
FeaturePlot(matrix.su, features = c("Flt1", "Anxa1", "Uba6", "Bnip3","Tceb2", "Plac1")) # EPC


message("+-------------------------------------------------------------------------------")
message("+                             Find Markers                                      ")
message("+-------------------------------------------------------------------------------")

matrix.su@meta.data$SCT_snn_res.1 <- factor(matrix.su@meta.data$SCT_snn_res.1, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18" , "19", "20", "21", "22", "23", "24", "25", "26" ))
Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.1


matrix.su.AllMarkers <- FindAllMarkers(matrix.su, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
#saveRDS(matrix.su.AllMarkers, "CTR_dscj1_MergedRealigned52_30L_SCT_matrix.su.AllMarkers.Rds")
matrix.su.AllMarkers <- readRDS("CTR_dscj1_MergedRealigned52_30L_SCT_matrix.su.AllMarkers.Rds")






message("--------------------------------------------------------------------------------")
message("+            Find Markers for all cluster pairs res 1.25                        ")
message("+-------------------------------------------------------------------------------")

resDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/MARKERS_SCT_res.1.25"
setwd(resDir)


head(matrix.su@meta.data,2)
Idents(matrix.su) <-  matrix.su@meta.data$SCT_snn_res.1.25
cluster_numbers <- c(0:(length(unique(matrix.su@meta.data$SCT_snn_res.1))-1))
marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = length(unique(matrix.su@meta.data$SCT_snn_res.1.25)), nrow = length(unique(matrix.su@meta.data$SCT_snn_res.1.25))))

# set up variables:
last_cluster <-length(unique(matrix.su@meta.data$SCT_snn_res.1.25))-1
l2fc_cutoff <- 0.25
reso <- "SCT_snn_res.1.25"
  
library(purrr)

findMarkers_for_cluster_pair <- function(i){
  if (j == i){
    print("same cluster")
  } else {
    tmp_markers <- FindMarkers(matrix.su, ident.1 = j, ident.2= i, logfc.threshold = l2fc_cutoff, verbose = TRUE)
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
saveRDS(result, paste0(reso, j, "_result_", "_l2fc" , l2fc_cutoff ,  ".rds"))
result_09 <- result

cluster_numbers <- c(10:last_cluster)
j <- "10"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_10 <- result

cluster_numbers <- c(11:last_cluster)
j <- "11"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_11 <- result

cluster_numbers <- c(12:last_cluster)
j <- "12"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_12 <- result

cluster_numbers <- c(13:last_cluster)
j <- "13"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_13 <- result

cluster_numbers <- c(14:last_cluster)
j <- "14"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_14 <- result

cluster_numbers <- c(15:last_cluster)
j <- "15"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_15 <- result

cluster_numbers <- c(16:last_cluster)
j <- "16"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_16 <- result

cluster_numbers <- c(17:last_cluster)
j <- "17"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_17 <- result

cluster_numbers <- c(18:last_cluster)
j <- "18"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_18 <- result

cluster_numbers <- c(19:last_cluster)
j <- "19"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_19 <- result

cluster_numbers <- c(20:last_cluster)
j <- "20"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_20 <- result

cluster_numbers <- c(21:last_cluster)
j <- "21"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_21 <- result

cluster_numbers <- c(22:last_cluster)
j <- "22"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_22 <- result

cluster_numbers <- c(23:last_cluster)
j <- "23"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_23 <- result

cluster_numbers <- c(24:last_cluster)
j <- "24"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_24 <- result

cluster_numbers <- c(25:last_cluster)
j <- "25"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_25 <- result

cluster_numbers <- c(26:last_cluster)
j <- "26"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_26 <- result

cluster_numbers <- c(27:last_cluster)
j <- "27"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_27 <- result

cluster_numbers <- c(28:last_cluster)
j <- "28"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_28 <- result

cluster_numbers <- c(29:last_cluster)
j <- "29"
result = map(cluster_numbers, safely_findMarkers)
saveRDS(result, paste0(reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
result_29 <- result


#savehistory(file = "SCT_marker_pairs_res1.0.Rhistory")




message("--------------------------------------------------------------------------------")
message("+                               RES 1.0 :::::                                  ")
message("+-------------------------------------------------------------------------------")

Project <- "CTR_dscj1_MergedRealigned52_30L_MARKERS"
resDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/MARKERS_SCT_res.1"
setwd(resDir)

last_cluster <- length(unique(matrix.su@meta.data$SCT_snn_res.1))-1
l2fc_cutoff <- 0.6
reso <- "SCT_snn_res.1"

marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = last_cluster+1, nrow = last_cluster+1 ))
rownames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20","c21","c22","c23","c24","c25","c26")
colnames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20","c21","c22","c23","c24","c25","c26")


result_00 <- readRDS("SCT_snn_res.1result_0_l2fc0.6.rds")
result_01 <- readRDS("SCT_snn_res.1result_1_l2fc0.6.rds") 
result_02 <- readRDS("SCT_snn_res.1result_2_l2fc0.6.rds") 
result_03 <- readRDS("SCT_snn_res.1result_3_l2fc0.6.rds") 
result_04 <- readRDS("SCT_snn_res.1result_4_l2fc0.6.rds") 
result_05 <- readRDS("SCT_snn_res.1result_5_l2fc0.6.rds") 
result_06 <- readRDS("SCT_snn_res.1result_6_l2fc0.6.rds") 
result_07 <- readRDS("SCT_snn_res.1result_7_l2fc0.6.rds") 
result_08 <- readRDS("SCT_snn_res.1result_8_l2fc0.6.rds") 
result_09 <- readRDS("SCT_snn_res.1result_9_l2fc0.6.rds") 
result_10 <- readRDS("SCT_snn_res.1result_10_l2fc0.6.rds") 
result_11 <- readRDS("SCT_snn_res.1result_11_l2fc0.6.rds") 
result_12 <- readRDS("SCT_snn_res.1result_12_l2fc0.6.rds") 
result_13 <- readRDS("SCT_snn_res.1result_13_l2fc0.6.rds") 
result_14 <- readRDS("SCT_snn_res.1result_14_l2fc0.6.rds") 
result_15 <- readRDS("SCT_snn_res.1result_15_l2fc0.6.rds") 
result_16 <- readRDS("SCT_snn_res.1result_16_l2fc0.6.rds") 
result_17 <- readRDS("SCT_snn_res.1result_17_l2fc0.6.rds") 
result_18 <- readRDS("SCT_snn_res.1result_18_l2fc0.6.rds") 
result_19 <- readRDS("SCT_snn_res.1result_19_l2fc0.6.rds") 
result_20 <- readRDS("SCT_snn_res.1result_20_l2fc0.6.rds") 
result_21 <- readRDS("SCT_snn_res.1result_21_l2fc0.6.rds") 
result_22 <- readRDS("SCT_snn_res.1result_22_l2fc0.6.rds") 
result_23 <- readRDS("SCT_snn_res.1result_23_l2fc0.6.rds") 
result_24 <- readRDS("SCT_snn_res.1result_24_l2fc0.6.rds") 
result_25 <- readRDS("SCT_snn_res.1result_25_l2fc0.6.rds") 
result_26 <- readRDS("SCT_snn_res.1result_26_l2fc0.6.rds") 

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
                    result_23 = result_23,
                    result_24 = result_24,
                    result_25 = result_25,
                    result_26 = result_26)
names(result_list)





message("--------------------------------------------------------------------------------")
message("+                               RES 1.25 :::::                                  ")
message("+-------------------------------------------------------------------------------")

Project <- "CTR_dscj1_MergedRealigned52_30L_MARKERS_res.1.25"
resDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/MARKERS_SCT_res.1.25"
setwd(resDir)

last_cluster <- length(unique(matrix.su@meta.data$SCT_snn_res.1.25))-1
l2fc_cutoff <- 0.6
reso <- "SCT_snn_res.1.25"

marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = last_cluster+1, nrow = last_cluster+1 ))
rownames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20","c21","c22","c23","c24","c25","c26", "c27", "c28", "c29")
colnames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20","c21","c22","c23","c24","c25","c26", "c27", "c28", "c29")


result_00 <- readRDS("SCT_snn_res.1.25_result_0_l2fc0.6.rds")
result_01 <- readRDS("SCT_snn_res.1.25_result_1_l2fc0.6.rds") 
result_02 <- readRDS("SCT_snn_res.1.25_result_2_l2fc0.6.rds") 
result_03 <- readRDS("SCT_snn_res.1.25_result_3_l2fc0.6.rds") 
result_04 <- readRDS("SCT_snn_res.1.25_result_4_l2fc0.6.rds") 
result_05 <- readRDS("SCT_snn_res.1.25_result_5_l2fc0.6.rds") 
result_06 <- readRDS("SCT_snn_res.1.25_result_6_l2fc0.6.rds") 
result_07 <- readRDS("SCT_snn_res.1.25_result_7_l2fc0.6.rds") 
result_08 <- readRDS("SCT_snn_res.1.25_result_8_l2fc0.6.rds") 
result_09 <- readRDS("SCT_snn_res.1.25_result_9_l2fc0.6.rds") 
result_10 <- readRDS("SCT_snn_res.1.25_result_10_l2fc0.6.rds") 
result_11 <- readRDS("SCT_snn_res.1.25_result_11_l2fc0.6.rds") 
result_12 <- readRDS("SCT_snn_res.1.25_result_12_l2fc0.6.rds") 
result_13 <- readRDS("SCT_snn_res.1.25_result_13_l2fc0.6.rds") 
result_14 <- readRDS("SCT_snn_res.1.25_result_14_l2fc0.6.rds") 
result_15 <- readRDS("SCT_snn_res.1.25_result_15_l2fc0.6.rds") 
result_16 <- readRDS("SCT_snn_res.1.25_result_16_l2fc0.6.rds") 
result_17 <- readRDS("SCT_snn_res.1.25_result_17_l2fc0.6.rds") 
result_18 <- readRDS("SCT_snn_res.1.25_result_18_l2fc0.6.rds") 
result_19 <- readRDS("SCT_snn_res.1.25_result_19_l2fc0.6.rds") 
result_20 <- readRDS("SCT_snn_res.1.25_result_20_l2fc0.6.rds") 
result_21 <- readRDS("SCT_snn_res.1.25_result_21_l2fc0.6.rds") 
result_22 <- readRDS("SCT_snn_res.1.25_result_22_l2fc0.6.rds") 
result_23 <- readRDS("SCT_snn_res.1.25_result_23_l2fc0.6.rds") 
result_24 <- readRDS("SCT_snn_res.1.25_result_24_l2fc0.6.rds") 
result_25 <- readRDS("SCT_snn_res.1.25_result_25_l2fc0.6.rds") 
result_26 <- readRDS("SCT_snn_res.1.25_result_26_l2fc0.6.rds") 
result_27 <- readRDS("SCT_snn_res.1.25_result_27_l2fc0.6.rds") 
result_28 <- readRDS("SCT_snn_res.1.25_result_28_l2fc0.6.rds") 
result_29 <- readRDS("SCT_snn_res.1.25_result_29_l2fc0.6.rds") 
result_30 <- readRDS("SCT_snn_res.1.25_result_30_l2fc0.6.rds") 
result_31 <- readRDS("SCT_snn_res.1.25_result_31_l2fc0.6.rds") 

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
                    result_23 = result_23,
                    result_24 = result_24,
                    result_25 = result_25,
                    result_26 = result_26,
                    result_27 = result_27,
                    result_28 = result_28,
                    result_29 = result_29)
names(result_list)

result_list <- readRDS("res.1.25_result-52_list.Rds")



message("--------------------------------------------------------------------------------")
message("+    Now caculate number of DE in pairwise clusters for choser RES              ")
message("+-------------------------------------------------------------------------------")

result_list <- readRDS("res.1.25_result-52_list.Rds")

resDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/MARKERS_SCT_res.1.25/CSV_tables"
setwd(resDir)

l2fc_cutoff <- 0.25

  
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
  marker_tbl[z,] <- c( rep(NA, 30-length(marker_tbl_list[[z]]) ), marker_tbl_list[[z]])
}
marker_tbl

length(marker_tbl_list[[1]])
length(marker_tbl_list[[2]])
length(marker_tbl_list[[3]])
length(marker_tbl_list[[4]])

max(marker_tbl, na.rm = T) # 79




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
      write.csv(tmp_df, paste0("CTR_dscj1_NewAlign52" , reso, "_l2fc", l2fc_cutoff,"_Markers_tbl__clusters_", name_y, ".vs.", name_x, ".csv"))
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

ht3 = Heatmap(as.matrix(marker_tbl),  col = f3, row_title = "", column_title = "Markers between cluster pairs (res 1.25,  abs(l2fc) > 0.25)", show_row_names = TRUE, heatmap_legend_param = list(title = "Number of markers", legend_height = unit(8, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left",  show_heatmap_legend = TRUE) #  width = unit(10, "cm"),
ht3


pdf(paste("Fig__Pairwise_MARKERS", Project, "ComplexHeatmap",  reso, "l2fc0.25", ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

#write.csv(marker_tbl, "CTR_dscj1_MergedRealigned52_30L_MARKERS___Pairwise_cluster_marker_tbl__res.1_l2f0.6.csv", quote = FALSE)


 
FeaturePlot(matrix.su, features = c("Lgals3", "Rhox9", "Phlda2", "Sparc","Sct", "Cdkn1c")) 
FeaturePlot(matrix.su, features = c("", "")) 






message("--------------------------------------------------------------------------------")
message("+                     GO  using enrichR                               ")
message("+-------------------------------------------------------------------------------")

library("enrichR")
enrichR_DB <- as.data.frame(listEnrichrDbs())

l2fc_cutoff <- 0.6
#l2fc_cutoff <- 0.25

markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/GO_MARKERS_SCT_res.1.25_l2fc0.6"
#markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/GO_MARKERS_SCT_res.1.25_l2fc0.25"
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
        write.csv(enrichR_RES, paste0("NewAlign52_SCT__Res.1.25__CCregr_", db, "_l2fc", l2fc_cutoff,"_", comparison_name, ".csv"))
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
l2fc_cutoff <- 0.6
#l2fc_cutoff <- 0.25
database <- "KEGG_2019_Mouse"
# GO_Cellular_Component_2018, GO_Biological_Process_2018, GO_Molecular_Function_2018  BioCarta_2016  WikiPathways_2019_Mouse KEGG_2019_Mouse

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

#GOBP_matrix_l2fc0.6  <- GO_matrix
#GOBP_matrix_l2fc0.25 <- GO_matrix
#GOCC_matrix_l2fc0.6  <- GO_matrix
#GOCC_matrix_l2fc0.25 <- GO_matrix
#GOMF_matrix_l2fc0.25 <- GO_matrix
#GOMF_matrix_l2fc0.6  <- GO_matrix
#BioCarta_matrix_l2fc0.25 <- GO_matrix
#BioCarta_matrix_l2fc0.6  <- GO_matrix
#WP_matrix_l2fc0.25 <- GO_matrix
#WP_matrix_l2fc0.6  <- GO_matrix
#Kegg_matrix_l2fc_0.6 <- GO_matrix
#Kegg_matrix_l2fc_0.25 <- GO_matrix

GO_res_list <- list(GOBP_matrix_l2fc0.6=GOBP_matrix_l2fc0.6,
                    GOBP_matrix_l2fc0.25=GOBP_matrix_l2fc0.25,
                    GOCC_matrix_l2fc0.6=GOCC_matrix_l2fc0.6,
                    GOCC_matrix_l2fc0.25=GOCC_matrix_l2fc0.25,
                    GOMF_matrix_l2fc0.6=GOMF_matrix_l2fc0.6,
                    GOMF_matrix_l2fc0.25=GOMF_matrix_l2fc0.25,
                    BioCarta_matrix_l2fc0.6=BioCarta_matrix_l2fc0.6,
                    BioCarta_matrix_l2fc0.25=BioCarta_matrix_l2fc0.25,
                    WP_matrix_l2fc0.6=WP_matrix_l2fc0.6,
                    WP_matrix_l2fc0.25=WP_matrix_l2fc0.25,
                    Kegg_matrix_l2fc_0.6=Kegg_matrix_l2fc_0.6,
                    Kegg_matrix_l2fc_0.25=Kegg_matrix_l2fc_0.25)
#saveRDS(GO_res_list, "NewAlign52_SCT_res0.2_GO_res_list.Rds")
names(GO_res_list)

for(i in seq_along(GO_res_list)){
  GO_res_list[[i]] <- GO_res_list[[i]][-c(1,2),]
}
head(GO_res_list[[1]],2)

for(i in seq_along(GO_res_list)){
  write.csv(GO_res_list[[i]], paste( "NewAlign52_SCT_res_1.25_", names(GO_res_list)[[i]], ".csv", sep = "_"),  quote = FALSE)
}




names(GO_res_list)[11] <-"Kegg_matrix_l2fc0.6" 
names(GO_res_list)[12] <-"Kegg_matrix_l2fc0.25"

GO_matrix_for_plotting <- GO_res_list[[1]]
GO_res_name <- names(GO_res_list)[[1]]

#for(i in seq_along(GO_res_list)){
for(i in c(7,8)){
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
  
  dim(GO_matrix2)
  f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
  ht1 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = T, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = F, cluster_rows = F , row_names_side ="left", width = unit(ncol(GO_matrix2), "cm"),) # width = unit(140, "cm"),
  ht1
  
  pdf(paste( Project, "ComplexHeatmap", "Fig__ALL_Pairwise_", GO_res_name, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix2)/2+5), height=nrow(GO_matrix2))
  par(bg=NA)
  draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
}










#l2fc_cutoff <- 0.6
#GO_matrix_for_plotting <- GOBP_matrix_l2fc0.6
#GO_res_name <- "GOBP_matrix_l2fc0.6"
#GO_matrix_for_plotting <- GOCC_matrix_l2fc0.6
#GO_res_name <- "GOCC_matrix_l2fc0.6"
#GO_matrix_for_plotting <- GOMF_matrix_l2fc0.6
#GO_res_name <- "GOMF_matrix_l2fc0.6"

l2fc_cutoff <- 0.25
#GO_matrix_for_plotting <- GOBP_matrix_l2fc0.25
#GO_res_name <- "GOBP_matrix_l2fc0.25"
#GO_matrix_for_plotting <- GOCC_matrix_l2fc0.25
#GO_res_name <- "GOCC_matrix_l2fc0.25"
GO_matrix_for_plotting <- GOMF_matrix_l2fc0.25
GO_res_name <- "GOMF_matrix_l2fc0.25"



GO_matrix_molten <- reshape2::melt(GO_matrix_for_plotting)
GO_matrix_molten <- na.omit(GO_matrix_molten)
GO_matrix_molten <- as.data.frame(table(GO_matrix_molten[,1]))
GO_matrix_molten <- GO_matrix_molten[order(GO_matrix_molten$Freq, decreasing = T),]
dim(GO_matrix_molten) # 1163 2


GO_matrix <- unique(GO_matrix_for_plotting)
GO_matrix2 <- GO_matrix
#GOBP_matrix2[is.na(GOBP_matrix2) == TRUE] <- 1
GO_matrix2[1:5,1:5]
rownames(GO_matrix2) <- GO_matrix2[,1]
GO_matrix2 <- GO_matrix2[,-1]

GO_matrix2_clust1 <- colnames(GO_matrix2)
GO_matrix2_clust1 <- gsub(".vs.*", "", GO_matrix2_clust1)
GO_matrix2_clust2 <- colnames(GO_matrix2)
GO_matrix2_clust2 <- gsub(".*.vs.", "", GO_matrix2_clust2)
GO_matrix2_clust2 <- gsub("X", "", GO_matrix2_clust2)

GO_matrix2[1:5,1:5]

f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
ht1 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = T, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = F, cluster_rows = F , row_names_side ="left",width = unit(140, "cm"),)
ht1

pdf(paste( Project, "ComplexHeatmap", "Fig__ALL_Pairwise_", GO_res_name,  "res.1.25", "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=70, height=40)
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()



GO_matrix2_1 <- as.data.frame(t(GO_matrix2))
GO_matrix2_1[1:15,1:3]
dim(GO_matrix2_1)


# All I and R  ::: clust 15 is T0
GO_matrix2_1a <- GO_matrix2_1[grepl("15.vs", rownames(GO_matrix2_1)),]
GO_matrix2_1b <- GO_matrix2_1[grepl("vs.15", rownames(GO_matrix2_1)),]
GO_matrix2_clust15 <- rbind(GO_matrix2_1a, GO_matrix2_1b)


GO_matrix2_clust15x <- GO_matrix2_clust15[,which(unlist(lapply(GO_matrix2_clust15, function(x)!all(is.na(x)))))]
GO_matrix2_clust15x[is.na(GO_matrix2_clust15x)] <- 1
dim(GO_matrix2_clust15x)

GO_matrix2_clust15x$clust <- rownames(GO_matrix2_clust15x)


# for clust 15:::
#GO_matrix2_clust15x$clust <- factor(GO_matrix2_clust15x$clust, levels = c("15.vs.21", "15.vs.28","08.vs.15","10.vs.15","06.vs.15", "03.vs.15", "05.vs.15", "02.vs.15","15.vs.23","09.vs.15","15.vs.16" ,  "15.vs.24", "15.vs.27",  "15.vs.29"))
   
  
GO_matrix2_clust15x <- GO_matrix2_clust15x[order(GO_matrix2_clust15x$clust),]
GO_matrix2_clust15x <- GO_matrix2_clust15x[,1:length(colnames(GO_matrix2_clust15x))-1]
GO_matrix2_clust15x <- as.data.frame(t(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( " \\(GO:.*" ,  "" , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "regulation" ,  "reg." , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "dependent" ,  "dep." , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "positive" ,  "+" , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "negative" ,  "-" , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "pathway" ,  "path." , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "chemical" ,  "chem." , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "polymerase" ,  "pol" , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "cotranslational" ,  "cotransl." , rownames(GO_matrix2_clust15x))
rownames(GO_matrix2_clust15x) <- gsub( "biosynthetic" ,  "biosynth." , rownames(GO_matrix2_clust15x))




f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
ht4 = Heatmap(as.matrix(GO_matrix2_clust15x),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE ,row_names_side ="left", width = unit(15, "cm"), row_dend_side = "right")
ht4


pdf(paste(Project, "ComplexHeatmap","Fig__Pairwise_GO_all___to_clust15_T0__", GO_res_name, "res.1.25", "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=20, height=15)
par(bg=NA)
draw(ht4, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

















message("--------------------------------------------------------------------------------")
message("+                            T0_c18 vs T0_c5 vs T0_c21                          ")
message("+-------------------------------------------------------------------------------")

matrix.su@meta.data$T0_clusts <- paste0(matrix.su@meta.data$Experiment, "_c",  matrix.su@meta.data$SCT_snn_res.1.25)
unique(matrix.su@meta.data$T0_clusts)

Idents(matrix.su) <- matrix.su@meta.data$T0_clusts

T0_c18_vs_T0_c5 <- FindMarkers(matrix.su, ident.1 = "0.0_c18", ident.2= "0.0_c5", verbose = TRUE)
T0_c18_vs_T0_c5 <- T0_c18_vs_T0_c5[T0_c18_vs_T0_c5$p_val_adj <0.05,]
T0_c18_vs_T0_c21 <- FindMarkers(matrix.su, ident.1 = "0.0_c18", ident.2= "0.0_c21", verbose = TRUE)
T0_c18_vs_T0_c21 <- T0_c18_vs_T0_c21[T0_c18_vs_T0_c21$p_val_adj <0.05,]
T0_c5_vs_T0_c21 <- FindMarkers(matrix.su, ident.1 = "0.0_c5", ident.2= "0.0_c21", verbose = TRUE)
T0_c5_vs_T0_c21 <- T0_c5_vs_T0_c21[T0_c5_vs_T0_c21$p_val_adj <0.05,]

FeaturePlot(matrix.su, features = c("Rpl18-ps1", "Ubb", "Rpl37","Eef1a1","Rpl41","Tmsb10","Npm1","Ybx1","Gm10263")) 

matrix.su_t0 <- subset(matrix.su, Experiment == "0.0")
head(matrix.su_t0@meta.data,2)
df <- as.data.frame(table(matrix.su_t0@meta.data[,c("T0_clusts", "orig.ident2")]))
df <- subset(df, df$T0_clusts %in% c("0.0_c18", "0.0_c21", "0.0_c5"))
df <- dcast(df, formula =  orig.ident2 ~ T0_clusts)

df2 <- as.data.frame(table(matrix.su@meta.data$T0_clusts))
df2$Experiment <- gsub( "_.*", "", df2$Var1)
df2$Cluster <- gsub( ".*_", "", df2$Var1)
df2 <- df2[,-1]
df2 <- df2[,c(2,3,1)]
df2 <- dcast(df2, formula = Cluster ~Experiment )

FeaturePlot(matrix.su, features = c("Saa3", "Ccl2", "Gata2", "Gata3", "Klf17"))

