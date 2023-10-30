# 
rm(list=ls())
gc()
options(bitmapType='cairo')


library(cowplot)
library(Seurat)
library(ggrepel)
library(ggplot2)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)



message("--------------------------------------------------------------------------------")
message("+       load sample table and alignment rates                                   ")
message("+-------------------------------------------------------------------------------")

#sampleTable <- read.csv("sampleTable2.csv")
uniq.align.rates <- read.csv("uniq.align.rates.csv", header = F)
colnames(uniq.align.rates) <- c("Batch","SLX", "SampleName", "AlignmentRate")
uniq.align.rates$Batch2 <- gsub( "CTR_DropSeq_", "Batch", uniq.align.rates$Batch)
uniq.align.rates <- subset(uniq.align.rates, uniq.align.rates$Batch2 %in% c( "Batch0004", "Batch0006", "Batch0007", "Batch0008" , "Batch0011", "Batch0014"))

uniq.align.rates$sampleLabels <- gsub( "SLX-", "", uniq.align.rates$SampleName)
uniq.align.rates$sampleLabels <- gsub( "\\.[ATXGC]{8}\\.", "_", uniq.align.rates$sampleLabels)
uniq.align.rates$sampleLabels <- paste(uniq.align.rates$Batch2 , uniq.align.rates$sampleLabels, sep= "_")
  #  Batch0011_15982_N716
  
#write.csv(sampleTable, "sampleTable.csv")

sampleTable <- read.csv("sampleTable.csv")
sampleTable$uniq.align.rate <- uniq.align.rates[match(sampleTable$sampleLabels, uniq.align.rates$sampleLabels),]$AlignmentRate



message("--------------------------------------------------------------------------------")
message("+                      alignments in samples                                    ")
message("+-------------------------------------------------------------------------------")

sampleTable_0.50 <- subset(sampleTable, sampleTable$uniq.align.rate > 50)
sampleTable_0.40 <- subset(sampleTable, sampleTable$uniq.align.rate > 40)

table(sampleTable_0.50$CellType)
table(sampleTable_0.50[,c(4,5)])
table(sampleTable_0.40[,c(4,5)])
table(sampleTable[,c(4,5)])

# same for cells not samples here::::




message("--------------------------------------------------------------------------------")
message("+       load realigned matrix.su (all samples)                                  ")
message("+-------------------------------------------------------------------------------")


#matrix.su <- readRDS("CTR_dscj1_Realigned_50bp_matrix.su_alignRate40_300g.Rds")
matrix.su <- readRDS("CTR_dscj1_Realigned_50bp_matrix.su_alignRate40_300g_added2unknownSamples.Rds")
head(matrix.su@meta.data,2)

matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0011_15987_N705",]$Experiment <- "I.4"

#Batch_A - 51bp (0006/0008)
#Batch_B – 99bp (0007/0011)
#Batch_C – 131bp (0004/0014) 
matrix.su@meta.data$NewBatch <- matrix.su@meta.data$Batch
matrix.su@meta.data$NewBatch <- ifelse( matrix.su@meta.data$NewBatch == "Batch0006" | matrix.su@meta.data$NewBatch == "Batch0008", "A00", matrix.su@meta.data$NewBatch)
matrix.su@meta.data$NewBatch <- ifelse( matrix.su@meta.data$NewBatch == "Batch0011" | matrix.su@meta.data$NewBatch == "Batch0007", "0B0", matrix.su@meta.data$NewBatch)
matrix.su@meta.data$NewBatch <- ifelse( matrix.su@meta.data$NewBatch == "Batch0014" | matrix.su@meta.data$NewBatch == "Batch0004", "00C", matrix.su@meta.data$NewBatch)


unique(matrix.su@meta.data$orig.ident)

table(matrix.su@meta.data[,c(6,7)])
table(matrix.su@meta.data[,6])

View(matrix.su@meta.data)

nrow(matrix.su@meta.data)             # 296897 cells (305095 cells with 2 samples added)
min(matrix.su@meta.data$alignRate)    # 40.17%
min(matrix.su@meta.data$nFeature_RNA) # 300


# subset alignment rate to 50% foe a round up number:::
# matrix.su <- subset(matrix.su, subset = alignRate >= 50% )


matrix.su@meta.data$Sample <- matrix.su@meta.data$orig.ident
head(matrix.su@meta.data,2)
unique(matrix.su@meta.data$Sample) # 163 (some removed because of alignment rate too low!)
unique(matrix.su@meta.data$Batch)  # 6 batches as expected


DefaultAssay( matrix.su ) <- "RNA"







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


sampleTable_sc2bulk <- sampleTable_0.40
#sampleTable_sc2bulk <- sampleTable[sampleTable$Batch == "Batch0007",]
#sampleTable_sc2bulk <- sampleTable[sampleTable$sampleName %in% c("15305_N701","9339_N710", "9616_N710", "9339_N728","9616_N728","9339_N716","9616_N716","9339_N721","9616_N721","9339_N719","9616_N719","9339_N715","9616_N715","9339_N719", "9616_N719", "9339_N707", "9616_N707"), ]


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
#saveRDS(rld, "CTR_dscj1_Dropseq_newAlignment50bp_sc-to-bulk__rld.Rds")



message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

#rld <- readRDS("CTR_dscj1_Dropseq_newAlignment50bp_sc-to-bulk__rld.Rds")
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

Supp_Fig_1_A <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(Batch) ))) + 
  geom_point(size = 5 , alpha = 0.5) + 
  xlab(pc1lab) + ylab(pc2lab) + 
  scale_colour_manual(name="Batch", values = c("green", "blue", "violet", "red", "orange", "yellow"))  +
  theme_classic() 
  #geom_text_repel(aes(label=Batch), col = "black") 
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


sampleTable_0.40 <- sampleTable_0.40[order(sampleTable_0.40$CellType),]
sampleTable_0.40 <- sampleTable_0.40[order(sampleTable_0.40$Age),]

norm_mat_sel2 <- norm_mat_sel[ , match( sampleTable_0.40$sampleLabels, colnames(norm_mat_sel) )]
colnames(norm_mat_sel2) <- paste(colnames(norm_mat_sel2), sampleTable_0.40$CellType, sep = "_")
norm_mat_sel2[1:5,1:5]

plot(hclust(as.dist(1-cor(norm_mat_sel))))
#plot(hclust(as.dist(1-cor(t(norm_mat_sel[,-1])))))

dist_mat <- as.dist(1-cor(norm_mat_sel2))
dist_mat <- as.dist(cor(norm_mat_sel2))
corr_mat <- cor(norm_mat_sel2)
min(corr_mat)

library(ggfortify)
#ggplot2::autoplot(dist_mat) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot2::autoplot(corr_mat) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#corr_mat_boot <- boot::corr(norm_mat_sel2)
#ggplot2::autoplot(corr_mat_boot) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#GGally::ggcorr(norm_mat_sel2)
#autoplot(as.dist(1-cor(t(norm_mat_sel)))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



hmisc_mat <- Hmisc::rcorr(norm_mat_sel2, type=c("pearson"))		#Matrix of Correlations and P-values
hmisc_mat$P
hmisc_mat$r

library(circlize)
library(ComplexHeatmap)

#f1 = colorRamp2( c(0, 0.5, 0.8, 1), c("white", "lightskyblue",  "deepskyblue3", "blue4"), space = "RGB") 
f1 = colorRamp2( c(0, 0.7, 0.85, 1), c("white", "lightskyblue",  "deepskyblue3", "blue4"), space = "RGB") 

ht1 = Heatmap(as.matrix(hmisc_mat$r),  col = f1, name = "sc-to-Bulk",  row_title = "", column_title = "sc-to-Bulk corr", show_row_names = TRUE, heatmap_legend_param = list(title = "Correlation", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = T, cluster_rows = T , row_names_side ="left")
ht1

pdf(paste("ComplexHeatmap",  "newAlignment50bp_CLUST", "Corr_sc_to_bulk_2unknown.pdf", sep="_"), onefile=FALSE, width=30, height=30) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

ht1 = Heatmap(as.matrix(hmisc_mat$r),  col = f1, name = "sc-to-Bulk",  row_title = "", column_title = "sc-to-Bulk corr", show_row_names = TRUE, heatmap_legend_param = list(title = "Correlation", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = F, cluster_rows = F , row_names_side ="left")
ht1

pdf(paste("ComplexHeatmap",  "newAlignment50bp", "Corr_sc_to_bulk_2unknown.pdf", sep="_"), onefile=FALSE, width=30, height=30) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


#ggplot2::autoplot(hmisc_mat$r) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_brewer( type = "seq", palette = f1, direction = 1)





message("+-------------------------------------------------------------------------------")
message("+                  plot heatmap for unknown samples                             ")
message("+-------------------------------------------------------------------------------")

all_samples <- as.data.frame(rownames(hmisc_mat$r))

# compare unnown samples with other ones
idxs <- c(grep("N727", all_samples$`rownames(hmisc_mat$r)`), grep("N715", all_samples$`rownames(hmisc_mat$r)`), grep("N704", all_samples$`rownames(hmisc_mat$r)`),   which(all_samples$`rownames(hmisc_mat$r)` == "Batch0006_7635_N702_R.1") )


selected_hmisc_mat_r <- hmisc_mat$r[idxs,idxs]
dim(selected_hmisc_mat_r)

ht1 = Heatmap(as.matrix(selected_hmisc_mat_r),  col = f1, name = "sc-to-Bulk",  row_title = "", column_title = "sc-to-Bulk corr", show_row_names = TRUE, heatmap_legend_param = list(title = "Correlation", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = F, cluster_rows = F , row_names_side ="left", width = unit(14, "cm"), height = unit(14, "cm"))
ht1

ht2 = Heatmap(as.matrix(selected_hmisc_mat_r),  col = f1, name = "sc-to-Bulk",  row_title = "", column_title = "sc-to-Bulk corr", show_row_names = TRUE, heatmap_legend_param = list(title = "Correlation", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = T, cluster_rows = T , row_names_side ="left", width = unit(14, "cm"))
ht2


message("+-------------------------------------------------------------------------------")
message("+          compare duplicates and replicates with r2 values                     ")
message("+-------------------------------------------------------------------------------")

# Duplicates of Set A

Dups_setA <- c("Batch0004_9617_N720", "Batch0004_9617_N721", "Batch0004_9617_N724", "Batch0004_9617_N726", "Batch0006_7633_N719","Batch0006_7635_N719", "Batch0007_7637_N724")
Dups_setA_expt <- paste0(Dups_setA, "_0.0" ) 
selected_hmisc_mat_r <- hmisc_mat$r[Dups_setA_expt,Dups_setA_expt]
dim(selected_hmisc_mat_r)


# Duplicates of Set B

Dups_setB <- c("Batch0004_9617_N722", "Batch0004_9617_N719", "Batch0004_9617_N718", "Batch0004_9617_N727","Batch0004_9617_N729", "Batch0006_7633_N704","Batch0006_7635_N704")
Dups_setB_expt <- paste0(Dups_setB, "_0.0" ) 
selected_hmisc_mat_r <- hmisc_mat$r[Dups_setB_expt,Dups_setB_expt]
dim(selected_hmisc_mat_r)


# Duplicates of Set 1A

Dups_set1A <- c("Batch0006_7633_N723", "Batch0006_7635_N723", "Batch0008_9616_N723", "Batch0008_9339_N723", "Batch0011_15982_N723","Batch0014_9339_N723", "Batch0014_9616_N723")
Dups_set1A_expt <- paste0(Dups_set1A, "_I.1" ) 

selected_hmisc_mat_r <- hmisc_mat$r[Dups_set1A_expt,Dups_set1A_expt]
dim(selected_hmisc_mat_r)




selected_hmisc_mat_r <- hmisc_mat$r[c(Dups_set1A_expt,Dups_setB_expt,Dups_setA_expt),c(Dups_set1A_expt,Dups_setB_expt,Dups_setA_expt)]
dim(selected_hmisc_mat_r)



f2 = colorRamp2( c(0, 0.80,0.9, 0.95, 1), c("white", "white","lightskyblue",  "deepskyblue3", "blue4"), space = "RGB") 
ht3 = Heatmap(as.matrix(selected_hmisc_mat_r),  col = f2, name = "sc-to-Bulk",  row_title = "", column_title = "sc-to-Bulk duplicates", show_row_names = TRUE, heatmap_legend_param = list(title = "Correlation", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = T, cluster_rows = T , row_names_side ="left")
ht3





message("+-------------------------------------------------------------------------------")
message("+                     plot frequency of R2                              ")
message("+-------------------------------------------------------------------------------")

df_r2_molten <- melt(dfx[])

plt <- ggplot(df_r2_molten, aes(x = value)) + geom_histogram(bins = 100) + theme_classic() + 
  geom_vline(xintercept=0.95,color = "red") +   
  geom_text(aes(x=0.94, label="r2 = 0.95", y=1000), colour="red", angle=90, text=element_text(size=14)) +
  geom_vline(xintercept=0.9,color = "blue") +   
  geom_text(aes(x=0.89, label="r2 = 0.90", y=1000), colour="blue", angle=90, text=element_text(size=14)) 

pdf(paste("Freq_of_r2", "_.pdf", sep=""), width=10,height=8) #, type = "cairo")
par(bg=NA)
plt
dev.off()






message("+-------------------------------------------------------------------------------")
message("+                          initial Seurat processing                            ")
message("+-------------------------------------------------------------------------------")

# store mitochondrial percentage in object meta data
matrix.su <- PercentageFeatureSet(matrix.su, pattern = "^mt-", col.name = "percent.mt")
max(matrix.su@meta.data$percent.mt)
#View(matrix.su@meta.data)
dim(matrix.su) #  19860 296897

# Visualize QC metrics as a violin plot
Idents(matrix.su) <- matrix.su@meta.data$Batch
#VlnPlot(matrix.su, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# subset MT 
matrix.su <- subset(matrix.su, subset = percent.mt < 5)
dim(matrix.su) #  19860 285970

# run sctransform
#matrix.su <- SCTransform(matrix.su, vars.to.regress = "percent.mt", verbose = T)

matrix.su <- NormalizeData(matrix.su)
matrix.su <- FindVariableFeatures(matrix.su, nfeatures = 1000)
matrix.su <- ScaleData(matrix.su, vars.to.regress = "percent.mt")

# These are now standard steps in the Seurat workflow for visualization and clustering
matrix.su <- RunPCA(matrix.su)

ElbowPlot(matrix.su)

matrix.su <- RunUMAP(matrix.su, dims = 1:30)
DimPlot(matrix.su, label = TRUE) + NoLegend()

matrix.su <- FindNeighbors(matrix.su, dims = 1:30)
matrix.su <- FindClusters(matrix.su)
DimPlot(matrix.su, label = TRUE) + NoLegend()

VizDimLoadings(matrix.su, dims = 1:2, reduction = "pca")
DimPlot(matrix.su, reduction = "pca")




#saveRDS(matrix.su, "CTR_dscj1_Dropseq_NewAlign50bp_initQC_scale1000mvg.Rds")
matrix.su <- readRDS("CTR_dscj1_Dropseq_NewAlign50bp_initQC_scale1000mvg.Rds")


VizDimLoadings(matrix.su, dims = 1:2, reduction = "pca")
print(matrix.su[["pca"]], dims = 1:5, nfeatures = 10)

Idents(matrix.su) <- matrix.su@meta.data$Batch
Idents(matrix.su) <- matrix.su@meta.data$Experiment

DimPlot(matrix.su, reduction = "pca")

DimPlot(matrix.su, reduction = "umap")


ElbowPlot(matrix.su)











message("+-------------------------------------------------------------------------------")
message("+                         Calculate Correlation Matrix                          ")
message("+-------------------------------------------------------------------------------")


lm_eqn <- function(df){
  x <- df[,1]
  y <- df[,2]
  m <- lm( x ~ y, df);
  r2 = format(summary(m)$r.squared, digits = 3)
  return(r2)
}


calcCorr <- function(dds, rld, sample1, sample2){
  
  scatter.log2    <- as.data.frame( log2( 1+counts(dds, normalized=TRUE)[, c(sample1,sample2)] ) )
  scatter.rlog    <- as.data.frame( assay(rld)[, c(sample1, sample2)] )
  scatter.log2.r2 <- lm_eqn(scatter.log2)
  scatter.rlog.r2 <- lm_eqn(scatter.rlog)
  lab.rlog.x      <- (max(scatter.rlog[,1]) / 10 * 5)
  lab.rlog.y      <- (max(scatter.rlog[,2]) / 10 * 9.25)
  
  lm(scatter.rlog)[1]

  message(paste0(sample1, ":::", sample2, ":::", scatter.rlog.r2, ":::", scatter.log2.r2))

}

calcCorr(dds, rld, 2, 9)







message("+-------------------------------------------------------------------------------")
message("+   Correlation Plots")
message("+-------------------------------------------------------------------------------")

plotScatters <- function(dds, rld, sample1, sample2, SPLIT){
  
  scatter.log2    <- as.data.frame( log2( 1+counts(dds, normalized=TRUE)[, c(sample1,sample2)] ) )
  scatter.rlog    <- as.data.frame( assay(rld)[, c(sample1, sample2)] )
  scatter.log2.r2 <- lm_eqn(scatter.log2)
  scatter.rlog.r2 <- lm_eqn(scatter.rlog)
  lab.rlog.x      <- (max(scatter.rlog[,1]) / 10 * 5)
  lab.rlog.y      <- (max(scatter.rlog[,2]) / 10 * 9.25)
  
  lm(scatter.rlog)[1]
  
  scatter.rlog_anno <- scatter.rlog
  scatter.rlog_anno$ratio <- scatter.rlog_anno[[1]] / scatter.rlog_anno[[2]] 
  scatter.rlog_anno$anno <- "not"
  scatter.rlog_anno$anno <- ifelse( scatter.rlog_anno[[1]] > scatter.rlog_anno[[2]] + 2 |  scatter.rlog_anno[[2]] > scatter.rlog_anno[[1]] + 2 , "anno", "not") 
  
  
  scatter.rlog_anno$external_gene_name <- rownames(scatter.rlog_anno)
  labeldata.ann <- scatter.rlog_anno[scatter.rlog_anno$anno == "anno",]
  
  message(paste0(sample1, ":::", sample2, ":::", scatter.rlog.r2, ":::", scatter.log2.r2))
  elementTextSize <- 12
  if(SPLIT == "Sample"){ 
    SplitXno <- 2} else {SplitXno <- 3}
  
  plot <- ggplot(data=scatter.rlog, aes(x=scatter.rlog[,1], y=scatter.rlog[,2])) + 
    geom_point(size=1, alpha=0.25) + geom_smooth(method="lm", se=FALSE) + 
    annotate("text", x=lab.rlog.x, y=lab.rlog.y, label=paste("r^2 = ", scatter.rlog.r2, sep=""), colour='blue', size=3 ) + 
    annotate("segment", x = min(scatter.rlog[,1]), xend = max(scatter.rlog[,1]), y = min(scatter.rlog[,2]), yend = max(scatter.rlog[,2]), colour = "red") +
    xlab(paste("Sample: ", as.character(s2c[sample1,SplitXno]), sep="")) + ylab(paste("Sample: ", as.character(s2c[sample2,SplitXno]), sep="")) + theme(text = element_text(size=elementTextSize)) +
    geom_abline(intercept = 2, slope = 1, colour = "grey") +
    geom_abline(intercept = -2, slope = 1, colour = "grey") #+
  # geom_label_repel(data=labeldata.ann,aes(x=labeldata.ann[,1], y=labeldata.ann[,2], label=external_gene_name),colour='purple',  size=3 )  
  
  #  point.padding = unit(0.25, "lines"),segment.size = 1, segment.color = 'purple',  nudge_x =1, nudge_y=-0.5
  #ggtitle(paste(Project, " rlog transformed counts \nSample ", as.character(sampleTable[sample1,1]), " vs ", "Sample", as.character(sampleTable[sample2,1]), sep=""))
  #print({labeldata.ann$external_gene_name[labeldata.ann$external_gene_name %in% top_outliersx$gene] })
}


# ------ ------  ------ ------  ------ ------  ------ ------  ------ ------ 
colnames(dds)
head(s2c,2)
# T.0 x 3 samples
p.sample.1 <- plotScatters(dds, rld, 2, 9, SPLIT)
p.sample.2 <- plotScatters(dds, rld, 2, 6, SPLIT)
p.sample.3 <- plotScatters(dds, rld, 6, 9, SPLIT)
p.sample.4 <- plotScatters(dds, rld, 20, 12, SPLIT)
p.sample.5 <- plotScatters(dds, rld, 20, 14, SPLIT)
p.sample.6 <- plotScatters(dds, rld, 14, 12, SPLIT)
p.sample.7 <- plotScatters(dds, rld, 20, 2, SPLIT)
p.sample.8 <- plotScatters(dds, rld, 14, 6, SPLIT)
p.sample.9 <- plotScatters(dds, rld, 12, 9, SPLIT)

#setwd(resDir)
pdf(paste("Realigned_QC", "_scatter_rlog_counts.pdf", sep=""), width=12,height=12) #, type = "cairo")
par(bg=NA)
plot_grid( p.sample.1, p.sample.2, p.sample.3,  p.sample.4,p.sample.5, p.sample.6,p.sample.7,p.sample.8,p.sample.9, ncol = 3, nrow = 3 )
dev.off()

# R.24 x4  samples
p.sample.1 <- plotScatters(dds, rld, 10, 3, SPLIT)
p.sample.2 <- plotScatters(dds, rld, 11, 5, SPLIT)
p.sample.3 <- plotScatters(dds, rld, 10, 11, SPLIT)
p.sample.4 <- plotScatters(dds, rld, 3, 5, SPLIT)

p.sample.5 <- plotScatters(dds, rld, 21, 16, SPLIT)
p.sample.6 <- plotScatters(dds, rld, 22, 18, SPLIT)
p.sample.7 <- plotScatters(dds, rld, 21, 22, SPLIT)
p.sample.8 <- plotScatters(dds, rld, 16, 18, SPLIT)

p.sample.9 <- plotScatters(dds, rld, 21, 10, SPLIT)
p.sample.10 <- plotScatters(dds, rld, 22, 11, SPLIT)
p.sample.11 <- plotScatters(dds, rld, 16, 3, SPLIT)
p.sample.12 <- plotScatters(dds, rld, 18, 5, SPLIT)

pdf(paste("Realigned_QC", "_R.24_scatter_rlog_counts.pdf", sep=""), width=16,height=12) #, type = "cairo")
par(bg=NA)
plot_grid( p.sample.1, p.sample.2, p.sample.3,  p.sample.4,p.sample.5, p.sample.6,p.sample.7,p.sample.8,p.sample.9,p.sample.10,p.sample.11,p.sample.12, ncol = 4, nrow = 3 )
dev.off()

# I.24: x4 samples
p.sample.1 <- plotScatters(dds, rld, 1, 4, SPLIT)
p.sample.2 <- plotScatters(dds, rld, 1, 7, SPLIT)
p.sample.3 <- plotScatters(dds, rld, 1, 8, SPLIT)
p.sample.4 <- plotScatters(dds, rld, 4, 7, SPLIT)

p.sample.5 <- plotScatters(dds, rld, 19, 17, SPLIT)
p.sample.6 <- plotScatters(dds, rld, 19, 15, SPLIT)
p.sample.7 <- plotScatters(dds, rld, 19, 13, SPLIT)
p.sample.8 <- plotScatters(dds, rld, 15, 17, SPLIT)

p.sample.9 <- plotScatters(dds, rld, 19, 1, SPLIT)
p.sample.10 <- plotScatters(dds, rld, 17, 4, SPLIT)
p.sample.11 <- plotScatters(dds, rld, 15, 7, SPLIT)
p.sample.12 <- plotScatters(dds, rld, 13, 8, SPLIT)


pdf(paste("Realigned_QC", "_I.24_scatter_rlog_counts.pdf", sep=""), width=16,height=12) #, type = "cairo")
par(bg=NA)
plot_grid( p.sample.1, p.sample.2, p.sample.3,  p.sample.4,p.sample.5, p.sample.6,p.sample.7,p.sample.8,p.sample.9,p.sample.10,p.sample.11,p.sample.12, ncol = 4, nrow = 3 )
dev.off()


#  15305-n701 vs othe 701 samples of batch 0014. t24 too, check a few pairs, raw counts to vst.
head(s2c[,c(1,3,4,7)],20)

#s2c <- s2c[,c(1,3,4,7)]
#s2c$N_name <- gsub( ".*_N", "N", s2c$sampleName)
#head(s2c[s2c$N_name == "N701",],20)
#which(s2c$sampleName ==  "9339_N701")
which(s2c$sampleLabels ==  "Batch0014_9339_N701") # 11
which(s2c$sampleLabels ==  "Batch0014_9616_N701") # 12
which(s2c$sampleLabels ==  "Batch0014_15305_N701") # 13
which(s2c$sampleLabels ==  "Batch0008_9616_N701") # 9
which(s2c$sampleLabels ==  "Batch0008_9339_N701") # 8


p.N701_1 <- plotScatters(dds, rld, sample1=11, sample2=12, SPLIT )
p.N701_2 <- plotScatters(dds, rld, 11, 13, SPLIT )
p.N701_3 <- plotScatters(dds, rld, 12, 13, SPLIT )
p.N701_4 <- plotScatters(dds, rld, 2, 5, SPLIT )
p.N701_5 <- plotScatters(dds, rld, 8, 9, SPLIT )
p.N701_6 <- plotScatters(dds, rld, 5, 13, SPLIT )
p.N701_7 <- plotScatters(dds, rld, 9, 12, SPLIT )
p.N701_8 <- plotScatters(dds, rld, 5, 8, SPLIT )
p.N701_9 <- plotScatters(dds, rld, 2, 11, SPLIT )
p.N701_10 <- plotScatters(dds, rld, 8, 11, SPLIT )
p.N701_11 <- plotScatters(dds, rld, 1, 8, SPLIT )
p.N701_12 <- plotScatters(dds, rld, 4, 9, SPLIT )


pdf(paste("T0_sample_N701", "_scatter_rlog_counts_lab_outliers.pdf", sep=""), width=15,height=5) #, type = "cairo")
par(bg=NA)
plot_grid( p.N701_1, p.N701_2, p.N701_3,  ncol = 3, nrow = 1 )
dev.off()


pdf(paste("T0_sample_N701_betweenBatches", "_scatter_rlog_counts_lab_outliers2.pdf", sep=""), width=20,height=15) #, type = "cairo")
par(bg=NA)
plot_grid( p.N701_1, p.N701_2, p.N701_3,  p.N701_4, p.N701_5, p.N701_6, p.N701_7, p.N701_9, p.N701_10, p.N701_8, p.N701_11, p.N701_12, ncol = 4, nrow = 3 )
dev.off()




p.1.2 <- plotScatters(dds, rld, 5, 2)
p.1.3 <- plotScatters(dds, rld, 5, 3)
p.1.4 <- plotScatters(dds, rld, 5, 4)
p.1.5 <- plotScatters(dds, rld, 5, 1)
p.2.3 <- plotScatters(dds, rld, 2, 4)
p.2.4 <- plotScatters(dds, rld, 2, 3)
p.2.5 <- plotScatters(dds, rld, 2, 1)
p.3.4 <- plotScatters(dds, rld, 3, 4)
p.3.5 <- plotScatters(dds, rld, 3, 1)
p.4.5 <- plotScatters(dds, rld, 4, 1)


title <- ggdraw() + draw_label(paste(Project, " rlog transformed counts between batch pairs", sep=""), fontface='bold')
theme_set(theme_cowplot(font_size=12)) 
px.1  <- plot_grid( p.1.2, p.1.3, p.1.4, p.1.5,   labels=c("A", "B", "C", "D" ), ncol = 4, nrow = 1 )
px.2  <- plot_grid( p.2.3, p.2.4, p.2.5,   labels=c( "E", "F", "G", "" ), ncol = 4, nrow = 1 )
px.3  <- plot_grid( p.3.4, p.3.5,   labels=c( "H", "I", "", "" ), ncol = 4, nrow = 1 )
px.4  <- plot_grid( p.4.5,   labels=c( "J", "", "", "" ), ncol = 4, nrow = 1 )

plot_grid( px.1, px.2, px.3, px.4, ncol = 1, nrow = 4 )

png(paste("T0_batch", "_scatter_rlog_counts_lab_outliers.png", sep=""), width=5000,height=5000, type = "cairo")
par(bg=NA)
plot_grid( px.1, px.2, px.3, px.4, ncol = 1, nrow = 4 )
dev.off()



#pdf(paste(Project, "_DESeq2_scatter_rlog_counts.pdf", sep=""), width=10,height=7, onefile=FALSE)
#par(bg=NA)
#plot_grid(title, px.1, ncol=1, rel_heights=c(0.1, 1))
#dev.off()

png(paste(Project, "_DESeq2_scatter_rlog_counts1.png", sep=""), width = 1500, height = 1000 )
par(bg=NA)
plot_grid(title, px.1, ncol=1, rel_heights=c(0.1, 1))
dev.off()




colnames(dds)

which(s2c$sampleLabels ==  "Batch0008_9339_N728") # 27
which(s2c$sampleLabels ==  "Batch0008_9616_N728") # 29
which(s2c$sampleLabels ==  "Batch0014_9339_N728") # 24
which(s2c$sampleLabels ==  "Batch0014_9616_N728") # 12

which(s2c$sampleLabels ==  "Batch0008_9616_N716") # 18
which(s2c$sampleLabels ==  "Batch0008_9339_N716") # 14
which(s2c$sampleLabels ==  "Batch0014_9616_N716") # 6
which(s2c$sampleLabels ==  "Batch0014_9339_N716") # 1

which(s2c$sampleLabels ==  "Batch0008_9616_N721") # 28
which(s2c$sampleLabels ==  "Batch0008_9339_N721") # 26
which(s2c$sampleLabels ==  "Batch0014_9616_N721") # 11
which(s2c$sampleLabels ==  "Batch0014_9339_N721") # 5

p.N728_1 <- plotScatters(dds, rld, 27, 29, SPLIT ) # within batch 0008
p.N728_2 <- plotScatters(dds, rld, 24, 12, SPLIT ) # within batch 0014
p.N728_3 <- plotScatters(dds, rld, 24, 27, SPLIT ) # 14 vs 8
p.N728_4 <- plotScatters(dds, rld, 12, 29, SPLIT ) # 14 vs 8

p.N716_1 <- plotScatters(dds, rld, 1, 6, SPLIT )   # within batch 0014
p.N716_2 <- plotScatters(dds, rld, 18, 14, SPLIT ) # within batch 0008
p.N716_3 <- plotScatters(dds, rld, 6, 18, SPLIT )  # 14 vs 8
p.N716_4 <- plotScatters(dds, rld, 1, 14, SPLIT )  # 14 vs 8

p.N721_1 <- plotScatters(dds, rld, 28, 26, SPLIT )
p.N721_2 <- plotScatters(dds, rld, 11, 5, SPLIT )
p.N721_3 <- plotScatters(dds, rld, 11, 28, SPLIT )
p.N721_4 <- plotScatters(dds, rld, 5, 26, SPLIT )


pdf(paste("N728_N716_N721_", "_scatter_rlog_counts_lab_outliers2.pdf", sep=""), width=20,height=15) #, type = "cairo")
par(bg=NA)
plot_grid( p.N728_1, p.N728_2, p.N728_3,  p.N728_4, p.N716_1, p.N716_2, p.N716_3, p.N716_4, p.N721_1, p.N721_2, p.N721_3, p.N721_4, ncol = 4, nrow = 3 )
dev.off()




colnames(dds)

# t.0 ::
which(s2c$sampleLabels ==  "Batch0004_9617_N720") # 44
which(s2c$sampleLabels ==  "Batch0004_9617_N718") # 46
which(s2c$sampleLabels ==  "Batch0004_9617_N719") # 45

which(s2c$sampleLabels ==  "Batch0008_9339_N728") # 14
which(s2c$sampleLabels ==  "Batch0008_9616_N728") # 19
which(s2c$sampleLabels ==  "Batch0014_9339_N728") # 53
which(s2c$sampleLabels ==  "Batch0014_9616_N728") # 59

which(s2c$sampleLabels ==  "Batch0011_15988_N704") # 6
which(s2c$sampleLabels ==  "Batch0007_7637_N724") # 20

p.T0_1 <- plotScatters(dds, rld, 44, 46, SPLIT ) # within batch 0004
p.T0_2 <- plotScatters(dds, rld, 44, 14, SPLIT ) # 4 vs 8
p.T0_3 <- plotScatters(dds, rld, 44, 53, SPLIT ) # 4 vs 14
p.T0_4 <- plotScatters(dds, rld, 44, 6, SPLIT ) # 4 vs 11
p.T0_5 <- plotScatters(dds, rld, 44, 20, SPLIT ) # 4 vs 7


# I.24 :::
which(s2c$sampleLabels ==  "Batch0006_7633_N710") # 35
which(s2c$sampleLabels ==  "Batch0006_7633_N706") # 38
which(s2c$sampleLabels ==  "Batch0008_9339_N728") # 14
which(s2c$sampleLabels ==  "Batch0008_9616_N710") # 53
which(s2c$sampleLabels ==  "Batch0007_7637_N728") # 25
which(s2c$sampleLabels ==  "Batch0011_15982_N728") # 3
which(s2c$sampleLabels ==  "Batch0014_9339_N728") # 53
which(s2c$sampleLabels ==  "Batch0014_9339_N710") # 50
which(s2c$sampleLabels ==  "") # 53


p.I24_1 <- plotScatters(dds, rld, 50, 53, SPLIT )   # within batch 0014
p.I24_2 <- plotScatters(dds, rld, 16, 14, SPLIT ) # within batch 0008
p.I24_3 <- plotScatters(dds, rld, 35, 38, SPLIT )  # within batch 0006
p.I24_5 <- plotScatters(dds, rld, 50, 16, SPLIT )   # 14 vs 8
p.I24_6 <- plotScatters(dds, rld, 50, 35, SPLIT ) #  14 vs 6
p.I24_7 <- plotScatters(dds, rld, 50, 25, SPLIT )  # 14 vs 7
p.I24_8 <- plotScatters(dds, rld, 50, 3, SPLIT )  # 14 vs 11
p.I24_9 <- plotScatters(dds, rld, 3, 25, SPLIT )  # 11 vs 7
p.I24_10 <- plotScatters(dds, rld, 3, 16, SPLIT )  # 11 vs 8
p.I24_11 <- plotScatters(dds, rld, 35, 16, SPLIT )  # 6 vs 8
p.I24_12 <- plotScatters(dds, rld, 25, 16, SPLIT )  # 7 vs 8
p.I24_13 <- plotScatters(dds, rld, 35, 3, SPLIT )  # 6 vs 11


pdf(paste("T0_I24_allBatches_", "_scatter_rlog_counts_lab_outliers2.pdf", sep=""), width=25,height=20) #, type = "cairo")
par(bg=NA)
plot_grid( p.T0_1,p.T0_2,p.T0_3,p.T0_4,p.T0_5, p.I24_1, p.I24_2, p.I24_3,p.I24_5, p.I24_6, p.I24_7, p.I24_8, p.I24_9, p.I24_10, p.I24_11, p.I24_12, p.I24_13, ncol = 5, nrow = 4 )
dev.off()









findOutliers <- function(dds, rld, sample1, sample2){
  
  scatter.log2    <- as.data.frame( log2( 1+counts(dds, normalized=TRUE)[, c(sample1,sample2)] ) )
  scatter.rlog    <- as.data.frame( assay(rld)[, c(sample1, sample2)] )
  scatter.log2.r2 <- lm_eqn(scatter.log2)
  scatter.rlog.r2 <- lm_eqn(scatter.rlog)
  lab.rlog.x      <- (max(scatter.rlog[,1]) / 10 * 5)
  lab.rlog.y      <- (max(scatter.rlog[,2]) / 10 * 9.25)
  
  lm(scatter.rlog)[1]
  
  scatter.rlog_anno <- scatter.rlog
  scatter.rlog_anno$ratio <- scatter.rlog_anno[[1]] / scatter.rlog_anno[[2]] 
  scatter.rlog_anno$anno <- "not"
  scatter.rlog_anno$anno <- ifelse( scatter.rlog_anno[[1]] > scatter.rlog_anno[[2]] + 2 |  scatter.rlog_anno[[2]] > scatter.rlog_anno[[1]] + 2 , "anno", "not") 
  
  
  scatter.rlog_anno$external_gene_name <- rownames(scatter.rlog_anno)
  labeldata.ann <- scatter.rlog_anno[scatter.rlog_anno$anno == "anno",]
  if(nrow(labeldata.ann) > 0 ){
    labeldata.ann$anno <- 1
    labeldata.ann$anno <- as.numeric(labeldata.ann$anno)
  }
  message(paste0(sample1, ":::", sample2, ":::", scatter.rlog.r2, ":::", scatter.log2.r2))
  
  colnames(labeldata.ann)[4] <- paste0( colnames(labeldata.ann)[1], "_",colnames(labeldata.ann)[2])
  outlier_df <- labeldata.ann[,c(4,5)]
  
}

o.1.2 <- findOutliers(dds, rld, 5, 2)
o.1.3 <- findOutliers(dds, rld, 5, 3)
o.1.4 <- findOutliers(dds, rld, 5, 4)
o.1.5 <- findOutliers(dds, rld, 5, 1)
o.2.3 <- findOutliers(dds, rld, 2, 4)
o.2.4 <- findOutliers(dds, rld, 2, 3)
o.2.5 <- findOutliers(dds, rld, 2, 1)
o.3.4 <- findOutliers(dds, rld, 3, 4)
o.3.5 <- findOutliers(dds, rld, 3, 1)
o.4.5 <- findOutliers(dds, rld, 4, 1)

outlier_df <- merge(o.1.2, o.1.3, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.1.4, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.1.5, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.2.3, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.2.4, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.2.5, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.3.4, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.3.5, by = "external_gene_name", all.x = TRUE, all.y = TRUE)
outlier_df <- merge(outlier_df, o.4.5, by = "external_gene_name", all.x = TRUE, all.y = TRUE)

rownames(outlier_df) <- outlier_df$external_gene_name
outlier_df_t <- as.data.frame(t(outlier_df[,-1]))

outlier_df <- outlier_df[,-1]
outlier_df[is.na(outlier_df)] <- 0

top_outliers <- as.data.frame(rowSums(outlier_df, na.rm = T))
colnames(top_outliers) <- "freq"
top_outliers$gene <- rownames(top_outliers)
top_outliers <- top_outliers[order(top_outliers$freq, decreasing = T),]
top_outliers <- subset(top_outliers, top_outliers$freq > 3)


for(i in seq_along(colnames(outlier_df))){
  outlier_df[[i]] <- as.numeric(as.character(outlier_df[[i]]))
}
outlier_genes <- as.data.frame(rowSums(outlier_df))
colnames(outlier_genes)[1] <- "freq"
outlier_genes$gene <- rownames(outlier_genes)
outlier_genes <- outlier_genes[order(outlier_genes$freq, decreasing = T),]

colnames(outlier_df)

outlier_0014 <- as.data.frame(rowSums(outlier_df[,c(1:4)]))
colnames(outlier_0014)[1] <- "freq"
outlier_0014$gene <- rownames(outlier_0014)
outlier_0014 <- outlier_0014[order(outlier_0014$freq, decreasing = T),]

outlier_0006 <- as.data.frame(rowSums(outlier_df[,c(4,7,9,10)]))
colnames(outlier_0006)[1] <- "freq"
outlier_0006$gene <- rownames(outlier_0006)
outlier_0006 <- outlier_0006[order(outlier_0006$freq, decreasing = T),]

outlier_0008 <- as.data.frame(rowSums(outlier_df[,c(2,6,8,9)]))
colnames(outlier_0008)[1] <- "freq"
outlier_0008$gene <- rownames(outlier_0008)
outlier_0008 <- outlier_0008[order(outlier_0008$freq, decreasing = T),]

outlier_0011 <- as.data.frame(rowSums(outlier_df[,c(3,5,8,10)]))
colnames(outlier_0011)[1] <- "freq"
outlier_0011$gene <- rownames(outlier_0011)
outlier_0011 <- outlier_0011[order(outlier_0011$freq, decreasing = T),]

outlier_0007 <- as.data.frame(rowSums(outlier_df[,c(1,5,6,7)]))
colnames(outlier_0007)[1] <- "freq"
outlier_0007$gene <- rownames(outlier_0007)
outlier_0007 <- outlier_0007[order(outlier_0007$freq, decreasing = T),]


