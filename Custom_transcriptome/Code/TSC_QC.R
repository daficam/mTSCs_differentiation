


library("dplyr")
library("methods")
library("utils")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Matrix")
library("matrixStats")
library("Seurat")
library("useful")
library("reshape2")
library("DESeq2")
library("biomaRt")
library("ComplexHeatmap")
library("circlize")




NUMCORES      <- 6
Project       <- "TSC_QC"


baseDir <- "/storage/CTR-Projects/CTR_DropSeq/TSC_Custom_Transcriptome/"
print(baseDir)


sampleFiles <- list.files(paste0(baseDir, "FeatureCounts"), 
                            pattern='*featureCounts_counts.txt$', recursive=T, full.names=T)

sampleNames <- gsub(".storage.*FeatureCounts.", "", sampleFiles)
sampleNames <- gsub("_R1_trimmed.bam_featureCounts_counts.txt", "", sampleNames)

sampleNames_Orig <- sampleNames
sampleNames_Orig <- gsub("$", "_R1", sampleNames_Orig)

sampleNames <- gsub("_L00[0-9]", "", sampleNames)

sampleTable <- data.frame(sampleNames=sampleNames, fileNameDGE=sampleFiles, origFiles=sampleNames_Orig)

# sampleTable$scComp <- sampleNames
# sampleTable$scComp <- gsub("lane5_PL_024_1",         "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane5_PL_024_2",         "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane5_PL_03_1",          "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane5_PL_03_2",          "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane5_PL_X24_1",         "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane5_PL_X24_2",         "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane5_PL_X3_1",          "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane5_PL_X3_2",          "SC", sampleTable$scComp)
# 
# sampleTable$scComp <- gsub("lane3_24h_DES_ACTTGAAT", "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane3_4d_DES_CGATGTAT",  "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane3_ctrl_CAGATCAT",    "SC", sampleTable$scComp)
# 
# sampleTable$scComp <- gsub("lane8_24-DES1_ACAGTGAT", "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane8_24-DES2_TAGCTTAT", "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane8_3_24_CGATGTAT",    "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane8_3-4D_GGCTACAT",    "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane8_3-ctrl_CTTGTAAT",  "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane8_4d-DES_GCCAATAT",  "SC", sampleTable$scComp)
# sampleTable$scComp <- gsub("lane8_ctrl_GATCAGAT",    "SC", sampleTable$scComp)
# 
# sampleTable$scComp <- gsub("lane.*",   "ALL", sampleTable$scComp)
# sampleTable$scComp <- gsub("Sample.*", "ALL", sampleTable$scComp)
# 
# 
# sampleTable <- subset(sampleTable, sampleTable$scComp=="SC")
# sampleTable <- sampleTable[grepl("PL",sampleTable$sampleNames),]
# 

print(sampleTable)
head(sampleTable)
nrow(sampleTable)

#write.csv(sampleTable[,c(1,3)], file=paste0(baseDir, "/", Project, "_BulkRNA_SampleList", ".csv"), quote=F, row.names=F)



ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'chromosome_name'), mart = ensembl)  
head(ensEMBL2id)
nrow(ensEMBL2id)


DESeqDataSetFromFeatureCounts <- function (sampleTable, directory = ".", design, ignoreRank = FALSE, ...) 
{
  # From https://www.biostars.org/p/277316/
  if (missing(design)) 
    stop("design is missing")
  l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory, fn), skip=2))
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1)))) 
    stop("Gene IDs (first column) differ between files.")
  tbl <- sapply(l, function(a) a$V7)
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[, 1]
  dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[, -(1:2), drop = FALSE], design = design, ignoreRank, ...)
  return(dds)
}






customPCA <- function(sampleTBL, RLD, TOPNUM, model, ensEMBL2id) {
  
  elementTextSize <- 12
  
  rv     <- rowVars(RLD)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sampleName, type=sampleTBL$scComp, pca$x)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, label=scores$sampleName, color=sampleTBL$scComp) ) +
                    geom_point(size = 3, alpha=0.5 ) + 
                    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=3, colour="black") +
                    scale_color_manual(name="Class", guide = 'none',
                                       values=c("SC"=rgb(0,0,255, maxColorValue=255),
                                                "ALL"=rgb(255,0,0, maxColorValue=255))) +
                    xlab(pc1lab) + ylab(pc2lab) + 
                    labs(title = paste0("Bulk RNA-Seq: Troph Stem Cells"),
                         subtitle = paste0("PCA Top ", TOPNUM, " MV"),
                         caption = "Data from Myriam Hemberger" )  +     
                    theme_bw() +
                    theme(text = element_text(size=elementTextSize), aspect.ratio=1, legend.position="none",
                          panel.background = element_blank(), axis.line = element_line(colour = "grey")) 
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=sampleTBL$scComp, label=scores$sampleName) ) +
    geom_point(size = 3, alpha=0.75 ) + 
    xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
  #  scale_shape_manual(name="Allele", values = c(17, 16, 15)) +
    labs(title = paste0("Bulk RNA-Seq: Troph Stem Cells"),
         subtitle = paste0("PCA Top ", TOPNUM, " MV"),
         caption = "Data from Myriam Hemberger" ) +
    theme(text = element_text(size=elementTextSize), legend.position="none") 
  
  loadings                 <- as.data.frame(pca$rotation)
  loadings$ensembl_gene_id <- rownames(loadings)
  loadings                 <- merge(loadings, ensEMBL2id, by="ensembl_gene_id")
  
  pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:50),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) + 
                          geom_point(size = 3 ) + xlab("") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:50),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + 
                          geom_point(size = 3 ) + xlab("") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca, plt.pca.nl, pca.1.25.plot, pca.2.25.plot) )
}


dds.Samples <- DESeqDataSetFromFeatureCounts(sampleTable=sampleTable, directory="", design= ~ 1   )
dds.Samples <- DESeq(dds.Samples, parallel=F)



message("+-------------------------------------------------------------------------------")
message("+ Table of Normalised Counts")
message("+-------------------------------------------------------------------------------")

normCounts                        <- counts(dds.Samples,      normalized=TRUE)
normCounts.df                     <- as.data.frame(normCounts)
normCounts.df$MeanNormCounts      <- rowMeans(normCounts.df) 
normCounts.df$Log2_MeanNormCounts <- log2(normCounts.df$MeanNormCounts+1)
normCounts.df$ensembl_gene_id     <- rownames(normCounts.df)
normCounts.df.annot               <- merge(normCounts.df,ensEMBL2id, by="ensembl_gene_id" )
normCounts.df.annot               <- normCounts.df.annot[ order(normCounts.df.annot$MeanNormCounts,decreasing=TRUE), ]
rownames(normCounts.df.annot)     <- normCounts.df.annot$ensembl_gene_id
normCounts.df.annot$description   <- gsub(" .Source.*", "", normCounts.df.annot$description)
normCounts.df.annot               <- normCounts.df.annot[, c(2:ncol(normCounts.df.annot))]
head(normCounts.df.annot)

#write.csv(normCounts.df.annot, file=paste0(baseDir, "/", Project, "_BulkRNA_normCounts", ".csv"), quote=F)



normCounts.df.annot.gender    <- normCounts.df.annot
normCounts.df.annot.gender    <- subset(normCounts.df.annot.gender, 
                                       (external_gene_name=="Xist" | external_gene_name=="Ddx3y" | 
                                        external_gene_name=="Kdm5d" | external_gene_name=="Uty" | 
                                        external_gene_name=="Zfy1"))

rownames(normCounts.df.annot.gender) <- normCounts.df.annot.gender$external_gene_name
normCounts.df.annot.gender           <- normCounts.df.annot.gender[ , !(names(normCounts.df.annot.gender) %in% 
                                                                   c("ensembl_gene_id","description","entrezgene","gene_biotype"))]
normCounts.df.annot.gender.mlt <- melt(normCounts.df.annot.gender)
head(normCounts.df.annot.gender.mlt)

plt.gender <- ggplot(normCounts.df.annot.gender.mlt, aes(x=external_gene_name, y=value, fill=chromosome_name)) + 
              geom_bar(stat="identity", position=position_dodge()) +
              facet_wrap(~ variable, ncol=4) +
              scale_fill_manual(name="Gender Linked Gene Sets", values=(c("X"="pink", "Y"="blue"))) +
              ylab("Normalised Read Counts") +
              xlab("Gender Linked Genes") +
              theme(text=element_text(size=8,  family="sans"),
              axis.text.y = element_text(size=8),
              axis.text.x = element_text(size=8, angle = 45, hjust = 1),
              legend.position="none")

#sampleTable$Gender <- c("F", "M", "F", "M", 
#                        "F", "F", "F", "F",
#                        "F", "F", "F", "F",
#                        "F", "M")


rld.Samples     <- vst(dds.Samples,    blind=F)

pca.plt.250     <- customPCA(sampleTable, assay(rld.Samples), 250, "rld.250",  ensEMBL2id)

pdf(paste0(baseDir, "/", Project, "_Fig.PCA.250_PC1PC2.pdf"),width=8,height=8)
par(bg=NA)
pca.plt.250[[1]]
dev.off()


pdf(paste0(baseDir, "/", Project, "_Fig.PCA.250.pdf"),width=10,height=7)
par(bg=NA)
plot_grid(pca.plt.250[[1]], pca.plt.250[[1]], plt.gender, pca.plt.250[[3]], nrow=2, ncol=2, scale = c(0.75, 1, 1, 0.75))
dev.off()

pc_ex           <- plot_grid(pca.plt.250[[3]], pca.plt.250[[4]], nrow=2)
pdf(paste0(baseDir, "/", Project, "_Fig.PCA.250.pdf"),width=14,height=7)
par(bg=NA)
plot_grid(pca.plt.250[[1]], pc_ex, ncol=2, labels=c("A","B"), rel_widths = c(1.25,1), align = "h", axis = "bt")
dev.off()

#
# Heatmap top 50
#

functionMakeHeatmaps <- function(MATRIX,TITLE) {
  
  base_mean = rowMeans(MATRIX)
  RNA.mat   = t(apply(MATRIX, 1, scale))
  colnames(RNA.mat) <- colnames(MATRIX)
  
  print( min(RNA.mat[,1]))
  print( max(RNA.mat[,2]))
  print(head(RNA.mat))

  col_RNA = colorRamp2(c(-2, 0, 2), c("darkgreen", "grey95", "purple"))
  hm.height <- nrow(RNA.mat)*.25
  
  ht2 = Heatmap(RNA.mat, name="RNA", col=col_RNA, cluster_columns=T, 
                show_row_names=T, show_column_names=T, width=unit(10,"cm"),
                row_names_gp = gpar(fontsize = 8), 
                column_title = paste0("TSC Bulk (PL Samples)\n", TITLE)) 
  ht.both <- ht2 
  
  return(ht.both)
}


normCounts.df.annot.var                 <- normCounts.df.annot
normCounts.df.annot.var$variance        <- rowVars( as.matrix(normCounts.df.annot.var[,c(1:8)]) )
normCounts.df.annot.var$ensembl_gene_id <- rownames(normCounts.df.annot.var)
normCounts.df.annot.var                 <- normCounts.df.annot.var[ order(normCounts.df.annot.var$variance,
                                                                          decreasing=T), ]
normCounts.df.annot.var                 <- normCounts.df.annot.var[c(1:50),]
rownames(normCounts.df.annot.var)       <- normCounts.df.annot.var$external_gene_name
normCounts.df.annot.var                 <- normCounts.df.annot.var[,c(1:8)]
normCounts.df.annot.var                 <- log2(normCounts.df.annot.var+1)
head(normCounts.df.annot.var)



hm <- functionMakeHeatmaps(normCounts.df.annot.var, "Top 50 Most Variable Genes")
pdf(paste0(baseDir, "/", Project, "_Fig.HeatMap_top50_mv.pdf"),width=6,height=10)
par(bg=NA)
hm 
dev.off()



normCounts.df.annot.sel                 <- normCounts.df.annot
normCounts.df.annot.sel                 <- subset(normCounts.df.annot.sel, 
                                                  (external_gene_name=="Cdx2"   | external_gene_name=="Eomes" | 
                                                   external_gene_name=="Esrrb"  | external_gene_name=="Elf5" | 
                                                   external_gene_name=="Gcm1"   | external_gene_name=="Syna" |  
                                                   external_gene_name=="Prl2c2" | external_gene_name=="Cx31" | 
                                                   external_gene_name=="Sox2"   | external_gene_name=="Bmp4" |   
                                                     external_gene_name=="Ovol2" |
                                                   external_gene_name=="Pcdh12" | external_gene_name=="Gjb3" ) )
rownames(normCounts.df.annot.sel)       <- normCounts.df.annot.sel$external_gene_name 
normCounts.df.annot.sel                 <- normCounts.df.annot.sel[,c(1:8)]
normCounts.df.annot.sel                 <- log2(normCounts.df.annot.sel+1)
print(normCounts.df.annot.sel)


hm2 <- functionMakeHeatmaps(normCounts.df.annot.sel, "Selected Marker Genes")
hm2

pdf(paste0(baseDir, "/", Project, "_Fig.HeatMap_markers.pdf"),width=6,height=7)
par(bg=NA)
hm2 
dev.off()


normCounts.df.annot.sel2                 <- normCounts.df.annot
normCounts.df.annot.sel2                 <- subset(normCounts.df.annot.sel2, 
                                                  (external_gene_name=="Cdx2"   | external_gene_name=="Eomes" | 
                                                     external_gene_name=="Esrrb"  | external_gene_name=="Elf5" | 
                                                     external_gene_name=="Fgfr2" |
                                                     external_gene_name=="Sox2"   | external_gene_name=="Nr0b1" ) )
rownames(normCounts.df.annot.sel2)       <- normCounts.df.annot.sel2$external_gene_name 
normCounts.df.annot.sel2                 <- normCounts.df.annot.sel2[,c(1:8)]
normCounts.df.annot.sel2                 <- log2(normCounts.df.annot.sel2+1)
print(normCounts.df.annot.sel2)


hm3 <- functionMakeHeatmaps(normCounts.df.annot.sel2, "Selected Marker Genes (Fig 1)")
hm3

pdf(paste0(baseDir, "/", Project, "_Fig.HeatMap_markers.Fig1.pdf"),width=6,height=7)
par(bg=NA)
hm2 
dev.off()



normCounts.df.annot.sel3                 <- normCounts.df.annot
normCounts.df.annot.sel3                 <- subset(normCounts.df.annot.sel3, 
                                                   (external_gene_name=="Cdx2"   | external_gene_name=="Eomes" | 
                                                      external_gene_name=="Esrrb"  | external_gene_name=="Elf5" | 
                                                      external_gene_name=="Nr0b1" | external_gene_name=="Spry4" |
                                                      external_gene_name=="Bmp4" | external_gene_name=="Mapk4" |
                                                      external_gene_name=="Dusp6" | 
                                                      external_gene_name=="Gcm1"   | external_gene_name=="Klf2" |
                                                      external_gene_name=="Sox2"   | external_gene_name=="Dppa1" ) )
rownames(normCounts.df.annot.sel3)       <- normCounts.df.annot.sel3$external_gene_name 
normCounts.df.annot.sel3                 <- normCounts.df.annot.sel3[,c(1:8)]
normCounts.df.annot.sel3                 <- log2(normCounts.df.annot.sel3+1)
print(normCounts.df.annot.sel3)

colnames((normCounts.df.annot.sel3))

hm3 <- functionMakeHeatmaps(normCounts.df.annot.sel3, "Selected Marker Genes (Fig 1d)")
hm3
