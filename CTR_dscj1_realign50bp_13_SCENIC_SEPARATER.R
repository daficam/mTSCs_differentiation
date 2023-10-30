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
})

options(width=200)
baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/"
setwd(baseDir)

#dataset <- "Harmony_500g_R"
dataset <- "Harmony_500g_I"
# harmony_T.24_500g
#dataset <- "harmony_T.24_500g"





message("--------------------------------------------------------------------------------")
message("+                        Initialize settings                                    ")
message("+-------------------------------------------------------------------------------")


if (dataset == "Harmony_500g_R"){
  Project <- "CTR_dscj1_HARMONY_500g_R_SCENIC"
  ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/SCENIC_R/"
  matrix.su_sep <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_R.Rds")
  setwd(ResDir)
  matrix.su_sep@meta.data$SCT_snn_res.0.8 <- factor(matrix.su_sep@meta.data$SCT_snn_res.0.8, levels= c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8"))
  Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.0.8
  cellInfo <- readRDS("cellInfo_metadata_harmony_500g_R_res0.8_ALL_CELLS.Rds")
  #exprMat <- GetAssayData(matrix.su_sep, assay = "SCT", slot = "data")
  #saveRDS(exprMat, "exprMat__harmony_R_500g.Rds")
  exprMat <- readRDS( "exprMat__harmony_R_500g.Rds")
} else if (dataset == "Harmony_500g_I"){
  Project <- "CTR_dscj1_HARMONY_500g_I_SCENIC"
  ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/SCENIC_I/"
  matrix.su_sep <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_500g_matrix.su_I.Rds")
  setwd(ResDir)
  matrix.su_sep@meta.data$SCT_snn_res.0.8 <- factor(matrix.su_sep@meta.data$SCT_snn_res.0.8, levels= c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0"))
  Idents(matrix.su_sep) <- matrix.su_sep@meta.data$SCT_snn_res.0.8
  #saveRDS(matrix.su_sep@meta.data, "cellInfo_metadata_harmony_500g_I_res0.8_ALL_CELLS.Rds")
  cellInfo <- readRDS("cellInfo_metadata_harmony_500g_I_res0.8_ALL_CELLS.Rds")
  #exprMat <- GetAssayData(matrix.su_sep, assay = "SCT", slot = "data")
  #saveRDS(exprMat, "exprMat__harmony_I_500g.Rds")
  exprMat <- readRDS( "exprMat__harmony_I_500g.Rds")
} else if (dataset == "harmony_T.24_500g"){
  Project <- "CTR_dscj1_HARMONY_500g_T24_SCENIC"
  ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g/SCENIC_T24"
  setwd(ResDir)
  #matrix.su_T24 <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g/CTR_dscj1_matrix.su_500g_T.24_SCT_HARMONY_clust.Rds")
  #Idents(matrix.su_T24) <- matrix.su_T24@meta.data$SCT_snn_res.0.4
  #Idents(matrix.su_T24) <- factor(Idents(matrix.su_T24), levels = c(1,2,3,9,7,8,0,4,5,6))
  #saveRDS(matrix.su_T24@meta.data, "cellInfo_metadata_harmony_500g_T24_res0.4_ALL_CELLS.Rds")
  cellInfo <- readRDS("cellInfo_metadata_harmony_500g_T24_res0.4_ALL_CELLS.Rds")
  #exprMat <- GetAssayData(matrix.su_T24, assay = "SCT", slot = "data")
  #saveRDS(exprMat, "exprMat__harmony_T.24_500g.Rds")
  exprMat <- readRDS( "exprMat__harmony_T.24_500g.Rds")
} else {"You didn't provide valid dataset!!"}



scenicOptions <- initializeScenic(datasetTitle= paste0("SCENIC", "_", Project), org="mgi", dbDir="cisTarget_databases", dbs=defaultDbNames[["mgi"]],  nCores=5)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 123

saveRDS(scenicOptions, file=paste0("int/scenicOptions_", Project ,".Rds") )





message("--------------------------------------------------------------------------------")
message("+                       Co-expression network                                   ")
message("+-------------------------------------------------------------------------------")

dim(exprMat)

genesKept <- geneFiltering(as.matrix(exprMat), scenicOptions)
length(genesKept) # 4663 for R out of 18396, 4738 for I.
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
saveRDS(exprMat_filtered, paste0("exprMat_filtered_", Project, ".Rds"))

exprMat_filtered <- as.matrix(exprMat_filtered)
exprMat_filtered[is.na(exprMat_filtered) ] <- 0

runCorrelation(as.matrix(exprMat_filtered), scenicOptions)


# Genie3 is very computationally extensive::
exprMat_log <- log2(exprMat_filtered+1)
# genie3 takes forever so do arboreto again!
#runGenie3(exprMat_log, scenicOptions)

exportsForArboreto(as.matrix(exprMat_log),scenicOptions, dir = "int")





message("--------------------------------------------------------------------------------")
message("+               Arboreto analysis in python                                     ")
message("+-------------------------------------------------------------------------------")

# now run Arboreto pythons script!!! after it gives output, transform it to Genie3 output and run next steps.
# need to produce: R_network_arboreto_numpy_output.tsv





message("--------------------------------------------------------------------------------")
message("+               Post-Arboreto analysis- step1                                   ")
message("+-------------------------------------------------------------------------------")
dataset <- "Harmony_500g_I"



if (dataset == "Harmony_500g_R"){
  Project <- "CTR_dscj1_HARMONY_500g_R_SCENIC"
  ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/SCENIC_R/"
  setwd(ResDir)
  scenicOptions <- readRDS("int/scenicOptions_CTR_dscj1_HARMONY_500g_R_SCENIC.Rds")
  cellInfo <- readRDS("cellInfo_metadata_harmony_500g_R_res0.8_ALL_CELLS.Rds")
  exprMat <- readRDS( "exprMat__harmony_R_500g.Rds")
  cluster_order <- c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8")
  #RGN_linklist <- read.table("Harmony_500g_R_network_arboreto_numpy_output.tsv", sep = "\t")
} else if (dataset == "Harmony_500g_I"){
  Project <- "CTR_dscj1_HARMONY_500g_I_SCENIC"
  ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/SCENIC_I/"
  setwd(ResDir)
  scenicOptions <- readRDS("int/scenicOptions_CTR_dscj1_HARMONY_500g_I_SCENIC.Rds")
  cellInfo <- readRDS("cellInfo_metadata_harmony_500g_I_res0.8_ALL_CELLS.Rds")
  exprMat <- readRDS( "exprMat__harmony_I_500g.Rds")
  #RGN_linklist <- read.table("Harmony_500g_I_network_arboreto_numpy_output.tsv", sep = "\t")
} else if (dataset == "harmony_T.24_500g"){
  Project <- "CTR_dscj1_HARMONY_500g_T24_SCENIC"
  ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g/SCENIC_T24"
  setwd(ResDir)
  scenicOptions <- readRDS("int/scenicOptions_CTR_dscj1_HARMONY_500g_T24_SCENIC.Rds")
  cellInfo <- readRDS("cellInfo_metadata_harmony_500g_T24_res0.4_ALL_CELLS.Rds")
  exprMat <- readRDS( "exprMat__harmony_T.24_500g.Rds")
  #RGN_linklist <- read.table("Harmony_500g_T24_network_arboreto_numpy_output.tsv", sep = "\t")
} else if (dataset == "harmony.orig_together_500g"){
  Project <- "CTR_dscj1_HARMONY_500g_Together_SCENIC"
  ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/SCENIC_together"
  setwd(ResDir)
  scenicOptions <- readRDS("int_SUBSETTED/scenicOptions_CTR_dscj1_HARMONY_orig_500g_SCENIC.Rds")
  cellInfo <- readRDS("SUBSETTED_FILES/cellInfo_metadata_64k_harmony.orig_500g_res0.8.Rds")
  exprMat <- readRDS( "SUBSETTED_FILES/exprMat_filtered_CTR_dscj1_HARMONY_orig_500g_SCENIC.Rds")
}


scenicOptions@inputDatasetInfo$int_01

dim(RGN_linklist)
head(RGN_linklist)
colnames(RGN_linklist) <- c("TF", "Target", "weight")
length(unique(RGN_linklist$TF))     #  R:  380
length(unique(RGN_linklist$Target)) #  R: 4663

corrMat <- readRDS("int/1.2_corrMat.Rds")
dim(corrMat)

RGN_linklist2 <- RGN_linklist[RGN_linklist$Target %in% rownames(corrMat),]
dim(RGN_linklist)
dim(RGN_linklist2)

saveRDS(RGN_linklist2, 'int/1.4_GENIE3_linkList.Rds')
#
runSCENIC_1_coexNetwork2modules(scenicOptions)
#08:53   Creating TF modules
#75%      90% 
#  1.275268 4.601185 
#Number of links between TFs and targets: 1574796
#[,1]
#nTFs          380
#nTargets     4663
#nGeneSets    2274
#nLinks    3469506


message("--------------------------------------------------------------------------------")
message("+               Post-Arboreto analysis- step2                                   ")
message("+-------------------------------------------------------------------------------")

nCores <- 12
library(BiocParallel) 
register(MulticoreParam(workers=nCores), default = TRUE) 
register(SnowParam(workers=nCores), default = TRUE)

scenicOptions@inputDatasetInfo$org <- "mgi" #"mgi"
scenicOptions@settings$db_mcVersion
scenicOptions@settings$dbs
motifAnnot <- getDbAnnotations(scenicOptions)

scenicOptions@settings$nCores <- 12
runSCENIC_2_createRegulons(scenicOptions, minGenes=10) #,  coexMethod=c("top10perTarget")) 

#12:04   Step 2. Identifying regulons
#tfModulesSummary:
  
#  top5perTarget top10perTarget          top50 top50perTarget           w001           w005 
#254            320            358            359            370            370 
#12:04   RcisTarget: Calculating AUC
#Scoring database:  [Source file: mm9-500bp-upstream-7species.mc9nr.feather]
#Scoring database:  [Source file: mm9-tss-centered-10kb-7species.mc9nr.feather]
#12:19   RcisTarget: Adding motif annotation
#Number of motifs in the initial enrichment: 1405144
#Number of motifs annotated to the matching TF: 15914
#12:22   RcisTarget: Prunning targets
#13:30   Number of motifs that support the regulons: 15914
#Preview of motif enrichment saved as: output/Step2_MotifEnrichment_preview.html
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0    12.0    67.0   221.7   370.0  1413.0 






message("--------------------------------------------------------------------------------")
message("+                step  3   -->   score cells                                    ")
message("+-------------------------------------------------------------------------------")

library("R2HTML")
library("doMC")


scenicOptions@settings$nCores <- 8 # not higher as otherwise will give error in score cells...
runSCENIC_3_scoreCells(scenicOptions, exprMat)

runSCENIC_4_aucell_binarize(scenicOptions)










message("--------------------------------------------------------------------------------")
message("+                   motif enrichment analysis                                   ")
message("+-------------------------------------------------------------------------------")

# top TSC/ diff markers ::: features = c("Elf5","Eomes","Cdx2","Gata3","Ascl2","Prdm1","Gcm1","Ovol2","Tfeb") )
# Trophoblast stem cell::: features = c("Elf5","Eomes","Cdx2","Sox2","Tead4","Esrrb") )
# Ectoplacental cone :::  features = c("Ascl2","Prdm1","Ets2") )
# Spongiotrophoblast/Junctional zone # Ascl2
# Trophoblast giant cells ::: features = c("Hand1","Stra13","Mdfi") , min.cutoff = 0)
# Labyrinthine trophoblast ::: features = c("Gcm1","Cebpa","Tfeb","Junb","Fosl1","Arnt", "Hif1a", "Hif1b", "Dlx3", "Pparg", "Pparbp") )

together_interesting_TFs <- c("Foxo3","Sp3","E2f8","Tead3","Grhl1","Ppard","Foxo4","Hdac2","Smarca4","Junb","Nfkb2","Creb3l2","Pou2f1","Mycn","Myc","Max","Hdac6","Srebf2","E2f7","Fos","Elf2","Bhlhe40", "Elf5", "Ddit3", "Irf7", "Gata1", "Atf4", "Stat2", "Msx2", "Mef2a", "Hif1a", "Mta3","Klf6", "Setdb1", "Jund", "Rest", "Ascl2", "Hes6","Atf6b", "Gata1", "Phf8", "Setdb1", "E2f8", "Gcm1", "Mef2a")

interesting_TFs <- unique(c("Elf5", "Sox2", "Ddit3" ,  "Foxo4", "Irf7", "Gata1",  "Atf4", "Bhlhe40",  "Stat2", "Gcm1", "Eomes", "Ets2", "Ascl2", "Hand1" , "Hif1a","Cdx2", "Gata3","Prdm1", "Ovol2", "Tfeb", "Tead4","Esrrb","Stra13", "Mdfi","Cebpa","Junb","Fosl1","Arnt","Hif1b", "Dlx3", "Pparg", "Pparbp" ))


#  EPC regulons : "Ascl2", "Hes6","Atf6b", "Gata1"
#  LB regulons:  "Phf8", "Setdb1", "E2f8", "Gcm1", "Mef2a"


interesting_TFs <- c(interesting_TFs, together_interesting_TFs)
# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
unique(motifEnrichment_selfMotifs_wGenes$highlightedTFs)

interesting_TFs[interesting_TFs %in% unique(motifEnrichment_selfMotifs_wGenes$highlightedTFs)]

tableSubset2 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs==interesting_TFs[[1]]]
viewMotifs(tableSubset2) 

#saveRDS(tableSubset, "I_tableSubset_Stat2.rds")
#saveRDS(motifEnrichment_selfMotifs_wGenes, "R_motifEnrichment_selfMotifs_wGenes.rds")




message("--------------------------------------------------------------------------------")
message("+                   regulon targets  analysis                                   ")
message("+-------------------------------------------------------------------------------")

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
#interesting_TFs <- c("Foxo3","Sp3","E2f8","Tead3","Grhl1","Ppard","Foxo4","Hdac2","Smarca4","Junb","Nfkb2","Creb3l2","Pou2f1","Mycn","Myc","Max","Hdac6","Srebf2","E2f7","Fos","Elf2","Bhlhe40", "Elf5", "Ddit3", "Foxo4", "Irf7", "Gata1", "Atf4", "Stat2")
tableSubset <- regulonTargetsInfo[TF %in% interesting_TFs,]
viewMotifs(tableSubset) 
head(regulonTargetsInfo)
head(tableSubset)



regulonTargetsInfo[gene %in% "Hif1a" & NES > 3]
regulonTargetsInfo[gene %in% "Hand1" & NES > 3]
regulonTargetsInfo[gene %in% "Stra13" & NES > 3 ]
regulonTargetsInfo[TF %in% "Elf5" & NES > 3]
regulonTargetsInfo[gene %in% "Eomes" & NES > 3]
regulonTargetsInfo[TF %in% "Gata3" & NES > 3 & highConfAnnot==TRUE]
regulonTargetsInfo[TF %in% "Sox2" & NES > 3]
regulonTargetsInfo[gene %in% "Ascl2" & NES > 3 ]
regulonTargetsInfo[TF %in% "Ets2" & NES > 3 ]
regulonTargetsInfo[TF %in% "Junb" & NES > 3 ]

regulonTargetsInfo[gene %in% "Gcm1"  ]

#write.csv(regulonTargetsInfo, paste0( Project, "_", "regulonTargetsInfo_ALL.csv"))









message("--------------------------------------------------------------------------------")
message("+                            plot AUC heatmap                                   ")
message("+-------------------------------------------------------------------------------")

### HEATMAP::: Regulators for clusters or known cell types

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUCx <- regulonAUC
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# for T24 reso is 0.4 not 0.8!!!
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$SCT_snn_res.0.8),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

head(regulonActivity_byCellType_Scaled)
#colnames(regulonActivity_byCellType_Scaled) <- c("0-R", "1-I", "2-I", "3-I", "4-R", "5-R", "6", "7", "8" ,"9-I")
rownames(regulonActivity_byCellType_Scaled)


all_TFs <- as.data.frame(rownames(regulonActivity_byCellType_Scaled))
colnames(all_TFs)[1] <- "TF"
all_TFs$TF_short <- all_TFs$TF
all_TFs$TF_short <- gsub( " \\(.*", "", all_TFs$TF_short)
all_TFs$TF_short <- gsub( "_extended", "", all_TFs$TF_short)

genes_for_heatmap <-  all_TFs$TF

regulonActivity_byCellType_Scaled2 <- regulonActivity_byCellType_Scaled[rownames(regulonActivity_byCellType_Scaled) %in% genes_for_heatmap,]

o = seriate(max(regulonActivity_byCellType_Scaled2) - regulonActivity_byCellType_Scaled2, method = "BEA_TSP")

set.seed(14)
ht1 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled2, #fontsize_row=3, 
                          color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                          treeheight_row=10, treeheight_col=10, border_color=NA, angle_col = "0", scale = "none", clustering_distance_cols= "correlation", row_order = get_order(o, 1), column_order = get_order(o, 2)) # 

ht1

pdf(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","SCT_snn_res.0.8_allTFRegulons","clustCorr.pdf", sep="_"), width=8,height=16)
par(bg=NA)
ht1
dev.off()

#regulonActivity_byCellType_Scaled2 <- regulonActivity_byCellType_Scaled2[,-which(colnames(regulonActivity_byCellType_Scaled2)== "20")]

ht2 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled2, #fontsize_row=3, 
                          color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                          treeheight_row=10, treeheight_col=10, border_color=NA, angle_col = "0", scale = "none", cluster_cols = FALSE)



pdf(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","SCT_snn_res.0.8_allTFRegulons","NotClust.pdf", sep="_"), width=8,height=16)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht2
dev.off()




# INHIBIT: 
EPC_clust <-   "3"
LP_clust  <-   "4"
T0_clust  <-   "9"

regulonActivity_3 <- regulonActivity_byCellType_Scaled2
colnames(regulonActivity_3)[colnames(regulonActivity_3)== EPC_clust] <- "EPC_clust" 
colnames(regulonActivity_3)[colnames(regulonActivity_3)== LP_clust] <- "LP_clust" 
colnames(regulonActivity_3)[colnames(regulonActivity_3)== T0_clust] <- "T0_clust" 
regulonActivity_3 <- as.data.frame(regulonActivity_3)
regulonActivity_3$ratio_EPC_LP <- regulonActivity_3$EPC_clust - regulonActivity_3$LP_clust

regulonActivity_3 <- regulonActivity_3[order(regulonActivity_3$ratio_EPC_LP),]
top_LP_regulons <- rownames(regulonActivity_3[c(1:10),])

regulonActivity_3 <- regulonActivity_3[order(regulonActivity_3$ratio_EPC_LP, decreasing = TRUE),]
top_EPC_regulons <- rownames(regulonActivity_3[c(1:10),])

regulonActivity_3 <- regulonActivity_3[order(regulonActivity_3$T0_clust, decreasing = TRUE),]
regulonActivity_3$maxValClust <- colnames(regulonActivity_3)[apply(regulonActivity_3,1,which.max)]
regulonActivity_3b <- regulonActivity_3[regulonActivity_3$maxValClust == "T0_clust",]
top_TSC_regulons <- rownames(regulonActivity_3b[ c(1:10),] )


regulonActivity_byCellType_Scaled3 <- regulonActivity_byCellType_Scaled2[rownames(regulonActivity_byCellType_Scaled2) %in% c(top_LP_regulons, top_EPC_regulons,top_TSC_regulons ),]


# monocle:
cluster_order_I <- c("9", "11" ,"16", "5" , "1" , "15", "6" , "12", "2" , "8" , "0" , "7" , "10", "14", "13", "4" , "3")
regulonActivity_byCellType_Scaled3 <- regulonActivity_byCellType_Scaled3[ , cluster_order_I]


ht2 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled3, #fontsize_row=3, 
                          color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                          treeheight_row=10, treeheight_col=10, border_color=NA, angle_col = "0", scale = "none", cluster_cols = TRUE)



pdf(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","top10_TSC_EPC_LP","Clust.pdf", sep="_"), width=6,height=6)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht2
dev.off()












#  EPC regulons : "Ascl2", "Hes6","Atf6b", "Gata1"
#  LB regulons:  "Phf8", "Setdb1", "E2f8", "Gcm1", "Mef2a"

pdf(paste( "UMAP_regulonAUC_SCENIC", dataset, "red_EomesElf5MycPou3f1__green_Ascl2Elf2Setdb1MycAtf6bGata1__blue_Gcm1E2f8Foxo4Hif1aMsx2Mef2aSrebf2Mta3RestPhf8",".pdf", sep="_"), width=15,height=11)
par(bg=NA)
regulonNames <- list(red=c("Eomes","Sox2", "Zfp553","Elf5","Pou3f1","Elf2"),
                     green=c( "Setdb1", "Myc","Ascl2","Atf6b","Gata1"),
                     blue=c("Foxo4", "", "Hif1a", "Msx2", "Mef2a","Srebf2","Mta3","Rest","Gcm1","E2f8","Phf8"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#text(-20, 5, attr(cellCol,"red"), col="red", cex=1.2, pos=4)
#text(-2,2, attr(cellCol,"blue"), col="blue", cex=1.2, pos=4)
#text(4,0, attr(cellCol,"green"), col="green3", cex=1.2, pos=4)
dev.off()


regulonNames <- c( "Elf2") # Pou3f1 Elf5
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)


message("--------------------------------------------------------------------------------")
message("+              replace tsne with seurat umap and binarize                       ")
message("+-------------------------------------------------------------------------------")


UMAPs <- as.data.frame(matrix.su_sep@reductions[["umap"]]@cell.embeddings)
#UMAPs$UMAP_1 <- UMAPs$UMAP_1 *-1
tsne_SCENIC <- readRDS("int/tSNE_AUC_50pcs_50perpl_from_SCENIC.Rds")
tsne_SCENIC$Y <- UMAPs
colnames(tsne_SCENIC$Y) <- c("tsne1" ,"tsne2")
#saveRDS(tsne_SCENIC, "int/tSNE_AUC_50pcs_50perpl.Rds")
#tsne_umap <- readRDS("int/tsne_to_umap.Rds")
#colnames(tsne_umap$Y) <- c("tsne1" ,"tsne2")
# on bash:::
#mv tSNE_AUC_50pcs_50perpl.Rds tSNE_AUC_50pcs_50perpl_orig_from_SCENIC.Rds
#mv tsne_to_umap.rds tSNE_AUC_50pcs_50perpl.Rds

scenicOptions@settings$defaultTsne

# setting up custom thresholds:::
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat, tSNE_fileName = "int/tSNE_AUC_50pcs_50perpl.Rds")
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
#saveRDS(newThresholds, "int/newThresholds.Rds" )
scenicOptions@fileNames$int


scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds" 
scenicOptions@settings$nCores <- 10 
# here rerunning using seurat umap instead of scenic tsne:::
runSCENIC_4_aucell_binarize(scenicOptions, skipTsne = TRUE, exprMat = exprMat)







message("--------------------------------------------------------------------------------")
message("+                     investigate and Plot regulons                             ")
message("+-------------------------------------------------------------------------------")
# together: Gcm1, Gata2, E2f8 !!, Mef2a?  
# Ascl2, Gata1? Atf6b?
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
regulon_names <- as.data.frame(cbind(onlyNonDuplicatedExtended(names(regulons))))
#regulon_names2 <- as.data.frame(cbind(names(regulons)))




pdf(paste( "UMAP_regulonAUC", dataset, "_green_E2f8", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( green=c("E2f8")), aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()

pdf(paste( "UMAP_regulonAUC", dataset, "_green_Elf3", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( green=c("Elf3")), aucType="AUC", aucMaxContrast=0.8,  cex=.4)
dev.off()


pdf(paste( "UMAP_regulonAUC", dataset, "_red_E2f8__blue_Elf3", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( red=c("E2f8"),blue=c("Elf3")), aucType="AUC", aucMaxContrast=0.6,  cex=.4)
dev.off()


pdf(paste( "UMAP_regulonAUC", dataset, "_green_Gcm1", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( green=c("Gcm1")), aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()


pdf(paste( "UMAP_regulonAUC", dataset, "_green_Ascl2", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( green=c("Ascl2")), aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()


regulonNames <- list(red=c("Eomes","Hand1"), #   "Elf5",      #"Cdx2", "Sox2", "Tead4", 
                     blue=c("Ascl2"), # , "Gata1" "Atf6b"  "Ets2", "Atf6b", "Gata1"  "Myc",
                     green=c("Gcm1"))  #  ,   "Phf8" "Mef2a", "Gata2", "Hif1a", "Srebf2", "Junb", "Hif1a", "Gata2"    # "Ovol2", "Pparg", 

pdf(paste( "UMAP_regulonAUC", dataset, "red_Eomes.Hand1__green_Gcm1__blue_Ascl2", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  regulonNames, aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()


regulonNames <- list(red=c("Eomes","Hand1"), #   "Elf5",      #"Cdx2", "Sox2", "Tead4", 
                     blue=c("Ascl2"))  #  ,   "Phf8" "Mef2a", "Gata2", "Hif1a", "Srebf2", "Junb", "Hif1a", "Gata2"    # "Ovol2", "Pparg", 

pdf(paste( "UMAP_regulonAUC", dataset, "red_Eomes.Hand1___blue_Ascl2", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  regulonNames, aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()


regulonNames <- list(red=c("Eomes", "Hand1", "Pou3f1"), #   "Elf5",      #"Cdx2", "Sox2", "Tead4", 
                     blue=c("Ascl2"), # , "Gata1" "Atf6b"  "Ets2", "Atf6b", "Gata1"  "Myc",
                     green=c("Gcm1", "E2f8"))  #  ,   "Phf8" "Mef2a", "Gata2", "Hif1a", "Srebf2", "Junb", "Hif1a", "Gata2"    # "Ovol2", "Pparg", 

pdf(paste( "UMAP_regulonAUC", dataset, "red_Eomes.Hand1.Pou3f1__green_E2f8.Gcm1__blue_Ascl2", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  regulonNames, aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()





pdf(paste( "UMAP_regulonAUC", dataset, "_green_Mef2a", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( green=c("Mef2a")), aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()


pdf(paste( "UMAP_regulonAUC", dataset, "_green_Phf8", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( green=c("Phf8")), aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()



pdf(paste( "UMAP_regulonAUC", dataset, "_blue_Gata2", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( blue=c("Gata2")), aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()

pdf(paste( "UMAP_regulonAUC", dataset, "_blue_Stat3", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( blue=c("Stat3")), aucType="AUC", aucMaxContrast=0.5,  cex=.4)
dev.off()
 
pdf(paste( "UMAP_regulonAUC", dataset, "_green_Atf6b", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( green=c("Atf6b")), aucType="AUC", aucMaxContrast=0.6,  cex=.4) # offColor = "gray50",
dev.off()

pdf(paste( "UMAP_regulonAUC", dataset, "_green_Myc", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
plotTsne_rgb(scenicOptions,  list( green=c("Myc")), aucType="AUC", aucMaxContrast=0.4,  cex=.4) 
dev.off()


#EPC lineage regulons : Ascl2, Myc, Atf6b, Gata1 (green), 
#LP lineage regulons : Gcm1, E2f8, Hif1a, Msx2, Mef2a, Srebf2, Mta3, Rest, Phf8 (blue) 
#TSC regulons: Eomes, Elf5, Pou3f1, Elf2 (red),


#regulons <- loadInt(scenicOptions, "regulons")
#names(regulons)
#regulons[c("Elf5_extended", "Elf3","Gata1", "Ets2", "Junb")] # "Gcm1","Esrrb"

regulonNames <- c( "Hand1", "E2f8", "Ascl2")
regulonNames <- c( "Elf5", "E2f8", "Ascl2")

regulonNames <- c( "Hand1", "Gcm1", "Ascl2")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(-18, 8, attr(cellCol,"red"), col="red", cex=1.2, pos=4)
text(8,12, attr(cellCol,"blue"), col="blue", cex=1.2, pos=4)
text(10,-9, attr(cellCol,"green"), col="green", cex=1.2, pos=4)


cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Elf5"), aucType="AUC", aucMaxContrast=0.8) #
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Eomes"), aucType="AUC", aucMaxContrast=0.8) #
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Hand1"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Taf1"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Kdm5a"), aucType="AUC", aucMaxContrast=0.8)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Pou3f1"), aucType="AUC", aucMaxContrast=0.8)



cellCol4 <- plotTsne_rgb(scenicOptions, regulonNames = c("Ets2"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Gcm1"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Irf7"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Stat2"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Gata2"), aucType="AUC", aucMaxContrast=0.8)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Gabpb1"), aucType="AUC", aucMaxContrast=0.8)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Phf8"), aucType="AUC", aucMaxContrast=0.8)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("E2f8"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Ddit3"), aucType="AUC", aucMaxContrast=0.6)


cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Hif1a"), aucType="AUC", aucMaxContrast=0.8)




# LP ::: Gcm1 Ddit3 Hif1a tbef==Atf3?/ Usf2 


cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Ascl2"), aucType="AUC", aucMaxContrast=0.6) #
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Gata1"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Atf6b"), aucType="AUC", aucMaxContrast=0.6)



cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Hif1a"), aucType="AUC", aucMaxContrast=0.6)

cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Bhlhe40"), aucType="AUC", aucMaxContrast=0.6)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Junb"), aucType="AUC", aucMaxContrast=0.8)


cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Usf2"), aucType="AUC", aucMaxContrast=0.8) # inhibits troph diff
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Atf3"), aucType="AUC", aucMaxContrast=0.8) # preeclamplsa gene
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Fos"), aucType="AUC", aucMaxContrast=0.8) # preeclamplsa gene


cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Irf7"), aucType="AUC", aucMaxContrast=0.6) # regulates Ly6e
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Stat2"), aucType="AUC", aucMaxContrast=0.6) # regulates Ly6e 

cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Irf1"), aucType="AUC", aucMaxContrast=0.6) # regulates Ly6e 

genes_LP <- c("Gcm1","Dlx3","Cebpa","Ovol2", "Junb", "Hif1a", "Tfeb") #  c("Gcm1","Cebpa","Tfeb","Junb","Fosl1","Arnt",  "Hif1b", "Dlx3", "Pparg", "Pparbp") 



cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Vezf1"), aucType="AUC")  
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Zfx"), aucType="AUC")  
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c(""), aucType="AUC")  
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c(""), aucType="AUC")  
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c(""), aucType="AUC")  
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c(""), aucType="AUC")  
cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c(""), aucType="AUC")  



#png(paste( "UMAP_regulonAUC", dataset, "_Eomes_Ascl2_Gcm1", "v1",".png", sep="_"), width=2000,height=2000, type= "cairo")
#par(bg=NA)
#cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#text(0, 10, attr(cellCol,"red"), col="red", cex=.7, pos=4)
#text(-20,-10, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
#dev.off()
Junb

genes_TSC <- c("Elf5","Esrrb","Cdx2","Eomes","Tead4","Sox2")
genes_GiantCells <- c("Mdfi","Prl8a9","Prl3d1","Prl3b1","Prl2c2","Ctsq")  # "Hand1","Stra13","Mdfi"
genes_JZ <- c("Pcdh12","Cdkn1c","Gjb3")
genes_EPC <- c("Ascl2", "Ets2") # "Hand1","Plet1","Prdm1",
genes_LP <- c("Gcm1","Dlx3","Cebpa","Ovol2", "Junb", "Hif1a", "Pparg") #  c("Gcm1","Cebpa","Tfeb","Junb","Fosl1","Arnt",  "Hif1b", "Dlx3", "Pparg", "Pparbp") 
genes_SynT <- c("Tfeb","Ly6e","Syna")
genes_Troph_ident <- c("Tfap2c","Gata3","Ets2","Hand1","Plet1")
genes_spiralArteryTGC <- "Prl2c2"
genes_canalTGC <- c("Prl3b1", "Prl2c2")
genes_parietalTGC <- c("Prl3b1","Prl2c2","Prl3d1")
interesting_genes <- c(genes_TSC, genes_GiantCells, genes_JZ, genes_EPC, genes_LP, genes_SynT, genes_Troph_ident, genes_spiralArteryTGC, genes_canalTGC, genes_parietalTGC)

interesting_genes[interesting_genes %in% c(regulonTargetsInfo$TF , regulonTargetsInfo$gene)]
#"Elf5"   "Eomes"  "Sox2"   "Cdkn1c" "Gjb3"   "Ascl2"  "Ets2"   "Junb"   "Hif1a"  "Ly6e"   "Gata3"   "Hand1" 
interesting_genes[interesting_genes %in% c(regulonTargetsInfo$TF )]
# "Elf5"  "Sox2"  "Ets2"  "Junb"  "Gata3"
unique(interesting_genes[interesting_genes %in% c(regulonTargetsInfo$gene )])
#  "Elf5"   "Eomes"  "Sox2"   "Cdkn1c" "Gjb3"   "Ascl2"  "Ets2"   "Junb"   "Hif1a"  "Ly6e"   "Gata3"  "Hand1" 


# top TSC/ diff markers ::: features = c("Elf5","Eomes","Cdx2","Gata3","Ascl2","Prdm1","Gcm1","Ovol2","Tfeb") )
# Trophoblast stem cell::: features = c("Elf5","Eomes","Cdx2","Sox2","Tead4","Esrrb") )
# Spongiotrophoblast/Junctional zone # Ascl2
# Labyrinthine trophoblast ::: features =)


# together: Gcm1, Gata2, E2f8 !!, Mef2a?
# Ascl2, Gata1? Atf6b?

regulonNames <- list(red=c("Eomes", "Cdx2", "Sox2", "Tead4", "Elf5", "Hand1"),
                     blue=c("Ascl2", "Ets2", "Atf6b", "Gata1"),
                     green=c("Gcm1", "Ovol2","Junb", "Hif1a", "Pparg", "E2f8", "Gata2"))
png(paste( "UMAP_regulonAUC", dataset, "_red.genes.TSC_blue.genes.EPC_Ascl2.Ets2.Gata1__green.genes_LP.Junb.Hif1a.E2f8.Gata2", "v4",".png", sep="_"), width=2000,height=1500, type= "cairo")
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC")
text(-20, 0, attr(cellCol,"red"), col="red", cex=1.5, pos=4)
text(-7,4, attr(cellCol,"blue"), col="blue", cex=1.5, pos=4)
text(-14,-3, attr(cellCol,"green"), col="green", cex=1.5, pos=4)
dev.off()


cellCol <- plotTsne_rgb(scenicOptions, regulonNames = c("Atf6b"), aucType="AUC", aucMaxContrast=0.8) # regulates Ly6e 

regulonNames <- list(red=c("Eomes", "Elf5", "Hand1"),
                     blue=c("Ascl2"), # "Gata1"
                     green=c("Gcm1", "E2f8")) # , "Gata2"
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC")
# for R:
text(-16, 8, attr(cellCol,"red"), col="red", cex=1.5, pos=4)
text(13, -5, attr(cellCol,"blue"), col="blue", cex=1.5, pos=4)
text(14,3, attr(cellCol,"green"), col="green", cex=1.5, pos=4)
# for I:
#text(-20, 7, attr(cellCol,"red"), col="red", cex=1.5, pos=4)
#text(7, 15, attr(cellCol,"blue"), col="blue", cex=1.5, pos=4)
#text(7,-10, attr(cellCol,"green"), col="green", cex=1.5, pos=4)



regulonNames <- list(red=c("Cdx2","Eomes","Tead4","Sox2", "Elf5"),
                     green=c("Ascl2", "Ets2"),
                     blue=c("Gcm1", "Ovol2", "Hif1a", "Junb"))
png(paste( "UMAP_regulonAUC", dataset, "red.genes_TSC_green.genes.Cdx2EomesTead4Sox2_EPC.Ascl2__blue.genes_LP.Gcm1Ovol2", "v4",".png", sep="_"), width=2000,height=1500, type= "cairo")
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.5)
text(-8, 0, attr(cellCol,"red"), col="black", cex=2.5, pos=4)
text(1,2, attr(cellCol,"green"), col="black", cex=2.5, pos=4)
text(4,3.5, attr(cellCol,"blue"), col="black", cex=2.5, pos=4)
dev.off()





regulonNames <- list( blue=c("Elf3"),
                     green=c("E2f8"))
png(paste( "UMAP_regulonAUC", dataset, "_red.genes.TSC_blue.genes.EPC_Ascl2.Ets2.Gata1__green.genes_LP.Junb.Hif1a.E2f8.Gata2", "v4",".png", sep="_"), width=2000,height=1500, type= "cairo")
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions,  list( blue=c("Elf3")), aucType="AUC")
#text(-20, 0, attr(cellCol,"red"), col="red", cex=1.5, pos=4)
#text(-7,4, attr(cellCol,"blue"), col="blue", cex=1.5, pos=4)
#text(-14,-3, attr(cellCol,"green"), col="green", cex=1.5, pos=4)
dev.off()







message("--------------------------------------------------------------------------------")
message("+                   Compare regulons between R and I                            ")
message("+-------------------------------------------------------------------------------")

# TEMPORA:::
cluster_order_I <- c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0")
cluster_order_R <- c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8")
# MONOCLE:::
cluster_order_I <- c("9", "11" ,"16", "5" , "1" , "15", "6" , "12", "2" , "8" , "0" , "7" , "10", "14", "13", "4" , "3")
cluster_order_R <- c("14", "13", "2" , "3" , "7", "0", "15", "6" , "12", "11", "5" , "10", "9" , "4" ,"8" ,"1")


setwd(baseDir)

regulonActivity_byCellType_Scaled_I <- as.data.frame(readRDS("SCENIC_R/regulonActivity_byCellType_Scaled___I.Rds"))
regulonActivity_byCellType_Scaled_R <- as.data.frame(readRDS("SCENIC_R/regulonActivity_byCellType_Scaled___R.Rds"))

Minval <- min(min(regulonActivity_byCellType_Scaled_I), min(regulonActivity_byCellType_Scaled_R))
Maxval <- max(max(regulonActivity_byCellType_Scaled_I), max(regulonActivity_byCellType_Scaled_R))

regulonActivity_byCellType_Scaled_I$regulon <- gsub( " .*", "", rownames(regulonActivity_byCellType_Scaled_I))
regulonActivity_byCellType_Scaled_R$regulon <- gsub( " .*", "", rownames(regulonActivity_byCellType_Scaled_R))

regulonActivity_byCellType_Scaled_I <- regulonActivity_byCellType_Scaled_I[order(regulonActivity_byCellType_Scaled_I$regulon),]
regulonActivity_byCellType_Scaled_R <- regulonActivity_byCellType_Scaled_R[order(regulonActivity_byCellType_Scaled_R$regulon),]

regulonActivity_byCellType_Scaled_I$regulon <- gsub( "_.*", "", regulonActivity_byCellType_Scaled_I$regulon)
regulonActivity_byCellType_Scaled_R$regulon <- gsub( "_.*", "", regulonActivity_byCellType_Scaled_R$regulon)

regulonActivity_byCellType_Scaled_I <- regulonActivity_byCellType_Scaled_I[!duplicated(regulonActivity_byCellType_Scaled_I$regulon),]
regulonActivity_byCellType_Scaled_R <- regulonActivity_byCellType_Scaled_R[!duplicated(regulonActivity_byCellType_Scaled_R$regulon),]

colnames(regulonActivity_byCellType_Scaled_I)[-ncol(regulonActivity_byCellType_Scaled_I)] <- paste0("c", colnames(regulonActivity_byCellType_Scaled_I)[-ncol(regulonActivity_byCellType_Scaled_I)])
#colnames(regulonActivity_byCellType_Scaled_R)[-ncol(regulonActivity_byCellType_Scaled_R)] <- paste0("c", colnames(regulonActivity_byCellType_Scaled_R)[-ncol(regulonActivity_byCellType_Scaled_R)])


regulonActivity_byCellType_Scaled_I <- regulonActivity_byCellType_Scaled_I[ , c( paste0("c", cluster_order_I),  "regulon")]
#regulonActivity_byCellType_Scaled_R <- regulonActivity_byCellType_Scaled_R[ , c( paste0("c", cluster_order_R),  "regulon")]

bothHTs <- merge(regulonActivity_byCellType_Scaled_I, regulonActivity_byCellType_Scaled_R, by = "regulon", all = TRUE)

rownames(bothHTs) <- bothHTs$regulon
bothHTs <- bothHTs[,-1]

regulonActivity_byCellType_Scaled_Ix <- bothHTs[,c(1:length(cluster_order_I))]
regulonActivity_byCellType_Scaled_Rx <- bothHTs[,c( (length(cluster_order_I)+1):ncol(bothHTs) )]

colnames(regulonActivity_byCellType_Scaled_Rx) <- paste0("c", colnames(regulonActivity_byCellType_Scaled_Rx))

library(seriation)
# find the right order as cannot cluster if missing values!
o = seriate(max(as.matrix(regulonActivity_byCellType_Scaled_I[,-18])) - as.matrix(regulonActivity_byCellType_Scaled_I[,-18]), method = "BEA_TSP")
get_order(o, 1) # row_order
get_order(o, 2) # col_order
tmp_I <- regulonActivity_byCellType_Scaled_I[get_order(o, 1),]
tmp_I$regulon # 107

o = seriate(max(as.matrix(regulonActivity_byCellType_Scaled_R[,-17])) - as.matrix(regulonActivity_byCellType_Scaled_R[,-17]), method = "BEA_TSP")
get_order(o, 1) # row_order
get_order(o, 2) # col_order
tmp_R <- regulonActivity_byCellType_Scaled_R[get_order(o, 1),]
tmp_R$regulon # 101

AllRegulons <- (c(regulonActivity_byCellType_Scaled_I$regulon, regulonActivity_byCellType_Scaled_R$regulon))
I_regulons <- tmp_I$regulon[!tmp_I$regulon %in% regulonActivity_byCellType_Scaled_R$regulon] # 20 unique for I
R_regulons <- tmp_R$regulon[!tmp_R$regulon %in% regulonActivity_byCellType_Scaled_I$regulon] # 14 unique for R
common_regulons <- tmp_I$regulon[tmp_I$regulon %in% regulonActivity_byCellType_Scaled_I$regulon & tmp_I$regulon %in% regulonActivity_byCellType_Scaled_R$regulon]

regulon_order <- c(common_regulons, rev(R_regulons) , I_regulons) 
split = c(c(rep("Common Regulons", length(common_regulons)), rep("Removal", length(R_regulons)), rep("Inhibit", length(I_regulons))))

regulonActivity_byCellType_Scaled_Ix <- regulonActivity_byCellType_Scaled_Ix[match( regulon_order, rownames(regulonActivity_byCellType_Scaled_Ix) ),]
regulonActivity_byCellType_Scaled_Rx <- regulonActivity_byCellType_Scaled_Rx[match( regulon_order, rownames(regulonActivity_byCellType_Scaled_Rx) ),]
regulonActivity_byCellType_Scaled_Rx <- regulonActivity_byCellType_Scaled_Rx[match( rownames(regulonActivity_byCellType_Scaled_Ix), rownames(regulonActivity_byCellType_Scaled_Rx) ),]


##### adding info if opposite regulons:::

col_mat_R <- rownames(regulonActivity_byCellType_Scaled_Rx)
col_mat_I <- rownames(regulonActivity_byCellType_Scaled_Ix)
col_mat_R <- ifelse(col_mat_R %in% c("Elf2","Nrf1","Gtf2f1","Myc","Kdm5b","Gata3","Junb"), "blue", "black")


#na_col = "lightgrey", column_names_gp = gpar(fontface = "italic", col = col_mat_R)
#na_col = "lightgrey", column_names_gp = gpar(fontface = "plain", col = col_mat_R)


f1 = colorRamp2( c(Minval,Minval/2,0,Maxval/2,Maxval), c("#0011FF", "#7F88FF",  "white", "#FF7F94", "#FF002A"), space = "LAB") 

ht1 = Heatmap(as.matrix(regulonActivity_byCellType_Scaled_Ix),  col = f1, name = "Inhibit",  row_title = "", column_title = "Inhibit regulons", show_row_names = TRUE, heatmap_legend_param = list(title = "AUC", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , width = unit(7, "cm"), row_names_side ="left", split = split,na_col = "lightgrey", row_names_gp = gpar(fontface = "plain", col = col_mat_R))

ht2 = Heatmap(as.matrix(regulonActivity_byCellType_Scaled_Rx),  col = f1, name = "Removal",  row_title = "", column_title = "Removal regulons", show_row_names = TRUE, heatmap_legend_param = list(title = "AUC", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(7, "cm"), row_names_side ="left", na_col = "lightgrey", row_names_gp = gpar(fontface = "plain", col = col_mat_R))

ht_list = ht1 + ht2
ht_list



pdf(paste("ComplexHeatmap_REGULONS_SPLIT_R_I", "colOrderMonocle_OppositeReg_blue.pdf", sep="_"), width=10,height=20)
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()




message("--------------------------------------------------------------------------------")
message("+   identify regulons changing in opposite direction between R and I            ")
message("+-------------------------------------------------------------------------------")


regulonActivity_byCellType_Scaled_Iy <- regulonActivity_byCellType_Scaled_Ix
regulonActivity_byCellType_Scaled_Ry <- regulonActivity_byCellType_Scaled_Rx

cluster_order_method <- "monocle"
#cluster_order_method <- "tempora"

if(cluster_order_method == "tempora"){
  # TEMPORA:::
  cluster_order_I <- c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0")
  cluster_order_R <- c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8")
  
} else {
  # MONOCLE:::
  cluster_order_I <- c("9", "11" ,"16", "5" , "1" , "15", "6" , "12", "2" , "8" , "0" , "7" , "10", "14", "13", "4" , "3")
  cluster_order_R <- c("14", "13", "2" , "3" , "7", "0", "15", "6" , "12", "11", "5" , "10", "9" , "4" ,"8" ,"1")
  
}



regulonActivity_byCellType_Scaled_Iy <- regulonActivity_byCellType_Scaled_Iy[, match( paste0( "c", cluster_order_I), colnames(regulonActivity_byCellType_Scaled_Iy) )]
regulonActivity_byCellType_Scaled_Ry <- regulonActivity_byCellType_Scaled_Ry[, match( paste0( "c", cluster_order_R), colnames(regulonActivity_byCellType_Scaled_Ry) )]




# mean regulon activity in top 5 diff clusters
regulonActivity_byCellType_Scaled_Iy$meanLate <- rowMeans(regulonActivity_byCellType_Scaled_Iy[,c(14:17) ])
regulonActivity_byCellType_Scaled_Ry$meanLate <- rowMeans(regulonActivity_byCellType_Scaled_Iy[,c(13:16) ])

regulonActivity_byCellType_Scaled_Iy$T0_clust <- regulonActivity_byCellType_Scaled_Iy[,1]
regulonActivity_byCellType_Scaled_Ry$T0_clust <- regulonActivity_byCellType_Scaled_Ry[,1]

#regulonActivity_byCellType_Scaled_Iy$T0_clust <- rowMeans(regulonActivity_byCellType_Scaled_Iy[,c(1:2) ])
#regulonActivity_byCellType_Scaled_Ry$T0_clust <- rowMeans(regulonActivity_byCellType_Scaled_Ry[,c(1:2) ])



regulonActivity_byCellType_Scaled_Iy$change <- regulonActivity_byCellType_Scaled_Iy$meanLate - regulonActivity_byCellType_Scaled_Iy$T0_clust
regulonActivity_byCellType_Scaled_Ry$change <- regulonActivity_byCellType_Scaled_Ry$meanLate - regulonActivity_byCellType_Scaled_Ry$T0_clust
regulonActivity_byCellType_Scaled_Iy <- regulonActivity_byCellType_Scaled_Iy[order(regulonActivity_byCellType_Scaled_Iy$change),]
regulonActivity_byCellType_Scaled_Ry <- regulonActivity_byCellType_Scaled_Ry[order(regulonActivity_byCellType_Scaled_Ry$change),]
regulonActivity_byCellType_Scaled_Iy$change_binary <- ifelse(regulonActivity_byCellType_Scaled_Iy$change > 0, "up_in_diff", "up_in_TSCs")
regulonActivity_byCellType_Scaled_Ry$change_binary <- ifelse(regulonActivity_byCellType_Scaled_Ry$change > 0, "up_in_diff", "up_in_TSCs")

regulonActivity__I <- regulonActivity_byCellType_Scaled_Iy[,c("T0_clust", "meanLate", "change","change_binary")]
regulonActivity__R <- regulonActivity_byCellType_Scaled_Ry[,c("T0_clust", "meanLate", "change","change_binary")]

colnames(regulonActivity__I) <- c("T0_clust_I", "meanLate_I", "change_I","change_binary_I")
colnames(regulonActivity__R) <- c("T0_clust_R", "meanLate_R", "change_R","change_binary_R")
binary_change_regulons <- merge(regulonActivity__I, regulonActivity__R, by = "row.names")
binary_change_regulons <- na.omit(binary_change_regulons)
binary_change_regulons_opposite <- binary_change_regulons[(binary_change_regulons$change_binary_I == "up_in_diff" & binary_change_regulons$change_binary_R == "up_in_TSCs") | (binary_change_regulons$change_binary_I == "up_in_TSCs" & binary_change_regulons$change_binary_R == "up_in_diff"),]

binary_change_regulons$diff_I_R <- binary_change_regulons$change_I - binary_change_regulons$change_R
#split <- binary_change_regulons_opposite

#binary_change_regulons_opposite <- binary_change_regulons[abs(binary_change_regulons$diff_I_R) > 1,]


f1 = colorRamp2( c(Minval,Minval/2,0,Maxval/2,Maxval), c("#0011FF", "#7F88FF",  "white", "#FF7F94", "#FF002A"), space = "LAB") 

mat_I <- regulonActivity_byCellType_Scaled_Iy[rownames(regulonActivity_byCellType_Scaled_Iy) %in% binary_change_regulons_opposite$Row.names, c(1:17)]

mat_R <- regulonActivity_byCellType_Scaled_Ry[rownames(regulonActivity_byCellType_Scaled_Ry) %in% binary_change_regulons_opposite$Row.names, c(1:16)]

mat_R <- mat_R[match( rownames(mat_I), rownames(mat_R) ),]


ht1 = Heatmap(mat_I,  col = f1, name = "Inhibit",  row_title = "", column_title = "Inhibit regulons", show_row_names = TRUE, heatmap_legend_param = list(title = "AUC", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(7, "cm"), row_names_side ="left")
#ht1

ht2 = Heatmap(mat_R,  col = f1, name = "Removal",  row_title = "", column_title = "Removal regulons", show_row_names = TRUE, heatmap_legend_param = list(title = "AUC", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(7, "cm"), row_names_side ="left")

ht_list = ht1 + ht2

ht_list


pdf(paste("ComplexHeatmap_REGULONS_oppositeChange", "colOrder", cluster_order_method, ".pdf", sep="_"), width=10,height=4)
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()




SELECTED_REGULONS <- c("Myc", "Nrf1", "Elf2", "Bmyc", "Bhlhe40", "Kdm5b", "Gtf2f1")


mat_I <- regulonActivity_byCellType_Scaled_Iy[rownames(regulonActivity_byCellType_Scaled_Iy) %in% SELECTED_REGULONS, c(1:17)]

mat_R <- regulonActivity_byCellType_Scaled_Ry[rownames(regulonActivity_byCellType_Scaled_Ry) %in% SELECTED_REGULONS, c(1:16)]

mat_R <- mat_R[match( rownames(mat_I), rownames(mat_R) ),]

ht1 = Heatmap(mat_I,  col = f1, name = "Inhibit",  row_title = "", column_title = "Inhibit regulons", show_row_names = TRUE, heatmap_legend_param = list(title = "AUC", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(7, "cm"), row_names_side ="left")
#ht1

ht2 = Heatmap(mat_R,  col = f1, name = "Removal",  row_title = "", column_title = "Removal regulons", show_row_names = TRUE, heatmap_legend_param = list(title = "AUC", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE , width = unit(7, "cm"), row_names_side ="left")

ht_list = ht1 + ht2
ht_list

pdf(paste("ComplexHeatmap_REGULONS_oppositeChange_SELECTED", "colOrderMonocle.pdf", sep="_"), width=10,height=3)
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()







message("--------------------------------------------------------------------------------")
message("+   check specific regulons           ")
message("+-------------------------------------------------------------------------------")

#regulonTargetsInfo_E2f8_TF <- regulonTargetsInfo[regulonTargetsInfo$TF =='Cebpg',]


regulonTargetsInfo_Phf8_TF <- regulonTargetsInfo[regulonTargetsInfo$TF =='Phf8',]


regulonTargetsInfo_Id2_Tfs <- regulonTargetsInfo[regulonTargetsInfo$gene == "Id2" & regulonTargetsInfo$highConfAnnot == "TRUE" & regulonTargetsInfo$nMotifs >=5 & regulonTargetsInfo$Genie3Weight >=1 , ]

regulonTargetsInfo_E2f8_TF <- regulonTargetsInfo[regulonTargetsInfo$TF =='E2f8',]

regulonTargetsInfo_E2f8_gene <- regulonTargetsInfo[ regulonTargetsInfo$gene == "E2f8",]


Klf6_targets_I <- regulonTargetsInfo[regulonTargetsInfo$TF == "Klf6",]
View(Klf6_targets_I)
Elf3_targets_I <- regulonTargetsInfo[regulonTargetsInfo$TF == "Elf3",]
Elf3_targets_I_conf <- regulonTargetsInfo[regulonTargetsInfo$TF == "Elf3" & regulonTargetsInfo$highConfAnnot == "TRUE",] # & regulonTargetsInfo$nMotifs >=5 & regulonTargetsInfo$Genie3Weight >=1 , ]
Elf3_targetedby_I <- regulonTargetsInfo[regulonTargetsInfo$gene == "Elf3",]
E2F8_targets_I <- regulonTargetsInfo[regulonTargetsInfo$TF == "E2f8",]
write.csv(Elf3_targets_I, "SCENIC_Inhibit_Elf3_targets.csv")
write.csv(Elf3_targets_I_conf, "SCENIC_Inhibit_Elf3_targets_HighConf.csv")


library("enrichR")
enrichR_DB <- as.data.frame(listEnrichrDbs())

database <- "WikiPathways_2019_Mouse"
# WikiPathways_2019_Mouse, GO_Biological_Process_2018, GO_Cellular_Component_2018, GO_Molecular_Function_2018,  BioCarta_2016, KEGG_2019_Mouse

SCENIC_targets <-  Elf3_targets_I_conf # regulonTargetsInfo_E2f8_TF #Klf6_targets_I

enrichR_RES <- as.data.frame(enrichr(SCENIC_targets$gene, databases = database))
enrichR_RES <- enrichR_RES[enrichR_RES$WikiPathways_2019_Mouse.Adjusted.P.value < 0.05,]
enrichR_RES <- enrichR_RES[,c(1,2,4,8,9)]

write.csv(enrichR_RES, "WikiPathways_SCENIC_Remove_Elf3_targets_HighConf.csv")
write.csv(enrichR_RES, "WikiPathways_SCENIC_Inhibit_Elf3_targets.csv")





