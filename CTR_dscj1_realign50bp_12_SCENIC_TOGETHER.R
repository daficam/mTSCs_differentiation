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
  
})

options(width=200)

Project <- "CTR_dscj1_HARMONY_orig_500g_SCENIC"


baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g"
setwd(baseDir)
ResDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/SCENIC_together"


#matrix.su <- readRDS("HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")
#sampleTable <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/sampleTable_MERGED_53.csv")

setwd(ResDir)



min(matrix.su@meta.data$nFeature_RNA)
max(matrix.su@meta.data$percent.mt) 
table(matrix.su@meta.data$SCT_snn_res.0.8)



message("--------------------------------------------------------------------------------")
message("+                                  Downsample                                   ")
message("+-------------------------------------------------------------------------------")
library(dplyr)

gc()
nrow(matrix.su@meta.data)

table(matrix.su@meta.data[,c("Experiment", "Batch")])
matrix.su@meta.data$Experiment_batch <- paste(matrix.su@meta.data$Experiment, matrix.su@meta.data$Batch, sep = "_")

tmp_tbl <- as.data.frame(matrix.su@meta.data$Experiment_batch)
rownames(tmp_tbl) <- rownames(matrix.su@meta.data)
colnames(tmp_tbl)[1] <- "Experiment"
tbl_order <- as.data.frame(table(tmp_tbl))
tbl_order <- tbl_order[order(tbl_order$Freq),]
Experiment_batch_names <-  as.character(tbl_order$tmp_tbl)

cells_to_keep <- list()
for (i in 1:length(Experiment_batch_names)) {
  sample_expt_tmp <- subset(tmp_tbl, tmp_tbl$Experiment == Experiment_batch_names[i])
  sample_expt_tmp$cells <- rownames(sample_expt_tmp)
  if (nrow(sample_expt_tmp) > 5000){
    cells_to_keepp <- sample_n(sample_expt_tmp,  5000 ) #  as.integer(nrow(sample_expt_tmp)*0.8)
    cells_to_keep[[i]] <- cells_to_keepp$cells
  } else {  cells_to_keep[[i]] <- sample_expt_tmp$cells  }
}
cells_to_keep <- unlist(cells_to_keep)
length(cells_to_keep)

matrix.su <- subset(matrix.su, cells = cells_to_keep)
matrix.su@meta.data$SCT_snn_res.0.8 <- factor(matrix.su@meta.data$SCT_snn_res.0.8, levels = c("12", "15", "0", "3", "17", "11", "13", "2", "14","20","18", "7", "8", "5", "9", "10", "16", "1", "19", "6", "4"))
Idents(matrix.su) <- matrix.su@meta.data$SCT_snn_res.0.8

#saveRDS(matrix.su, "matrix.su.downsample_64k_harmony.orig_500g_res0.8.Rds")
#cellInfo <- saveRDS(matrix.su@meta.data, "cellInfo_metadata_64k_harmony.orig_500g_res0.8.Rds")

saveRDS(matrix.su@meta.data, "cellInfo_metadata_ALL_harmony.orig_500g_res0.8.Rds")


#cellInfo <- saveRDS(matrix.su@meta.data, "cellInfo_metadata_ALL_harmony.orig_500g_res0.8.Rds")
#exprMat <- GetAssayData(matrix.su, assay = "SCT", slot = "data")
#saveRDS(exprMat, "exprMat_ALL_harmony.orig_500g_res0.8.Rds")


message("--------------------------------------------------------------------------------")
message("+                        Initialize settings                                    ")
message("+-------------------------------------------------------------------------------")


exprMat <- GetAssayData(readRDS("SUBSETTED_FILES/matrix.su.downsample_64k_harmony.orig_500g_res0.8.Rds"))
cellInfo <- readRDS("SUBSETTED_FILES/cellInfo_metadata_64k_harmony.orig_500g_res0.8.Rds")
#cellInfo <- readRDS("cellInfo_metadata_ALL_harmony.orig_500g_res0.8.Rds")

scenicOptions <- initializeScenic(datasetTitle= paste0("SCENIC", "_ALL_", Project), org="mgi", dbDir="cisTarget_databases", dbs=defaultDbNames[["mgi"]],  nCores=5)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 123

saveRDS(scenicOptions, file=paste0("scenicOptions_ALL_", Project ,".Rds") )




message("--------------------------------------------------------------------------------")
message("+                       Co-expression network                                   ")
message("+-------------------------------------------------------------------------------")

dim(exprMat)

# these below is when using ALL cells as too memory heavy:
#exprMat1 <- exprMat[,c(1:60000)]
#exprMat2 <- exprMat[,c(60001:ncol(exprMat))]
#genesKept1 <- geneFiltering(as.matrix(exprMat1), scenicOptions)
#genesKept2 <- geneFiltering(as.matrix(exprMat2), scenicOptions)
#genesKept <- c(genesKept1, genesKept2)
#exprMat_filtered <- exprMat[genesKept, ]


genesKept <- geneFiltering(as.matrix(exprMat), scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

saveRDS(exprMat_filtered, paste0("exprMat_filtered_ALL_", Project, ".Rds"))
#exprMat_filtered <- readRDS("exprMat_filtered_ALL_CTR_dscj1_HARMONY_orig_500g_SCENIC.Rds")


exprMat_filtered <- as.matrix(exprMat_filtered)
exprMat_filtered[is.na(exprMat_filtered) ] <- 0

runCorrelation(as.matrix(exprMat_filtered), scenicOptions)


# Genie3 is very computationally extensive::
exprMat_log <- log2(exprMat_filtered+1)
# genie3 takes forever so do arboreto again!
#runGenie3(exprMat_log, scenicOptions)

exportsForArboreto(as.matrix(exprMat_log), scenicOptions, dir = "int")





message("--------------------------------------------------------------------------------")
message("+               Arboreto analysis in python                                     ")
message("+-------------------------------------------------------------------------------")

# now run Arboreto pythons script!!! after it gives output, transform it to Genie3 output and run next steps.
# need to produce: R_network_arboreto_numpy_output.tsv




message("--------------------------------------------------------------------------------")
message("+               Post-Arboreto analysis- step1,2                                 ")
message("+-------------------------------------------------------------------------------")

setwd(ResDir)

cellInfo <- readRDS("SUBSETTED_FILES/cellInfo_metadata_64k_harmony.orig_500g_res0.8.Rds")
scenicOptions <- readRDS("int/scenicOptions_CTR_dscj1_HARMONY_orig_500g_SCENIC.Rds")



scenicOptions@inputDatasetInfo$int_01

RGN_linklist <- read.table("int/Harmony.orig.500g_network_arboreto_numpy_output.tsv", sep = "\t")
dim(RGN_linklist)
head(RGN_linklist)
colnames(RGN_linklist) <- c("TF", "Target", "weight")
length(unique(RGN_linklist$TF)) #[1] 378
length(unique(RGN_linklist$Target)) #[1]  4699

corrMat <- readRDS("int/1.2_corrMat.Rds")
dim(corrMat)

RGN_linklist2 <- RGN_linklist[RGN_linklist$Target %in% rownames(corrMat),]
dim(RGN_linklist)
dim(RGN_linklist2)
#saveRDS(RGN_linklist2, 'int/1.4_GENIE3_linkList.Rds')

runSCENIC_1_coexNetwork2modules(scenicOptions)



nCores <- 14
library(BiocParallel) 
register(MulticoreParam(workers=nCores), default = TRUE) 
register(SnowParam(workers=nCores), default = TRUE)

scenicOptions@inputDatasetInfo$org <- "mgi" #"mgi"
scenicOptions@settings$db_mcVersion
scenicOptions@settings$dbs
motifAnnot <- getDbAnnotations(scenicOptions)

scenicOptions@settings$nCores <- 14
runSCENIC_2_createRegulons(scenicOptions, minGenes=10,  coexMethod=c("top10perTarget")) 
# specified coexMethod not all as too memory intensive
#16:19   Step 2. Identifying regulons
#tfModulesSummary:
#  top10perTarget 
#320 
#16:19   RcisTarget: Calculating AUC
#Scoring database:  [Source file: mm9-500bp-upstream-7species.mc9nr.feather]
#Scoring database:  [Source file: mm9-tss-centered-10kb-7species.mc9nr.feather]
#16:21   RcisTarget: Adding motif annotation
#Number of motifs in the initial enrichment: 225398
#Number of motifs annotated to the matching TF: 3068
#16:21   RcisTarget: Prunning targets
#16:30   Number of motifs that support the regulons: 3068
#Preview of motif enrichment saved as: output/Step2_MotifEnrichment_preview.html
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    8.00   22.50   56.63   74.75  417.00 





message("--------------------------------------------------------------------------------")
message("+                step  3   -->   score cells                                    ")
message("+-------------------------------------------------------------------------------")

library("R2HTML")
library("doMC")


exprMat <- GetAssayData(readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/SCENIC_together/SUBSETTED_FILES/matrix.su.downsample_64k_harmony.orig_500g_res0.8.Rds"))

scenicOptions@settings$nCores <- 8 # not higher as otherwise will give error in score cells...
runSCENIC_3_scoreCells(scenicOptions, exprMat)

#runSCENIC_4_aucell_binarize(scenicOptions)




message("--------------------------------------------------------------------------------")
message("+              replace tsne with seurat umap and binarize                       ")
message("+-------------------------------------------------------------------------------")

matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/SCENIC_together/SUBSETTED_FILES/matrix.su.downsample_64k_harmony.orig_500g_res0.8.Rds")

UMAPs <- as.data.frame(matrix.su@reductions[["umap"]]@cell.embeddings)
UMAPs$UMAP_1 <- UMAPs$UMAP_1 *-1
tsne_SCENIC <- readRDS("int/tSNE_AUC_50pcs_50perpl_orig_from_SCENIC.Rds")
tsne_SCENIC$Y <- UMAPs
#saveRDS(tsne_SCENIC, "int/tsne_to_umap.Rds")
#tsne_umap <- readRDS("int/tsne_to_umap.Rds")
#colnames(tsne_umap$Y) <- c("tsne1" ,"tsne2")
# on bash:::
#mv tSNE_AUC_50pcs_50perpl.Rds tSNE_AUC_50pcs_50perpl_orig_from_SCENIC.Rds
#mv tsne_to_umap.rds tSNE_AUC_50pcs_50perpl.Rds

scenicOptions@settings$defaultTsne

# setting up custom thresholds:::
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat, tSNE_fileName = "int/tsne_to_umap.Rds")
savedSelections <- shiny::runApp(aucellApp)
newThresholds <- savedSelections$thresholds
#saveRDS(newThresholds, "int/newThresholds.Rds" )
scenicOptions@fileNames$int


scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds" 
scenicOptions@settings$nCores <- 10 
# here rerunning using seurat umap instead of scenic tsne:::
runSCENIC_4_aucell_binarize(scenicOptions, skipTsne = TRUE, exprMat = exprMat)





message("--------------------------------------------------------------------------------")
message("+                   motif enrichment analysis                                   ")
message("+-------------------------------------------------------------------------------")

interesting_TFs <- c("Elf5", "Sox2", "Ddit3" ,  "Foxo4", "Irf7", "Gata1",  "Atf4", "Bhlhe40",  "Stat2", "Klf6"  )

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
unique(motifEnrichment_selfMotifs_wGenes$highlightedTFs)
tableSubset2 <- motifEnrichment_selfMotifs_wGenes[highlightedTFs==interesting_TFs[[1]]]
viewMotifs(tableSubset2) 

#saveRDS(tableSubset, "I_tableSubset_Stat2.rds")
#saveRDS(motifEnrichment_selfMotifs_wGenes, "I_motifEnrichment_selfMotifs_wGenes.rds")



message("--------------------------------------------------------------------------------")
message("+                   regulon targets  analysis                                   ")
message("+-------------------------------------------------------------------------------")

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
interesting_TFs <- c("Foxo3","Sp3","E2f8","Tead3","Grhl1","Ppard","Foxo4","Hdac2","Smarca4","Junb","Nfkb2","Creb3l2","Pou2f1","Mycn","Myc","Max","Hdac6","Srebf2","E2f7","Fos","Elf2","Bhlhe40", "Elf5", "Ddit3", "Irf7", "Gata1", "Atf4", "Stat2", "Msx2", "Mef2a", "Hif1a", "Mta3","Klf6", "Setdb1", "Jund", "Rest")
tableSubset <- regulonTargetsInfo[TF %in% interesting_TFs,]
viewMotifs(tableSubset) 
head(regulonTargetsInfo)
head(tableSubset)


# top TSC/ diff markers ::: features = c("Elf5","Eomes","Cdx2","Gata3","Ascl2","Prdm1","Gcm1","Ovol2","Tfeb") )
# Trophoblast stem cell::: features = c("Elf5","Eomes","Cdx2","Sox2","Tead4","Esrrb") )
# Ectoplacental cone :::  features = c("Ascl2","Prdm1","Ets2") )
# Spongiotrophoblast/Junctional zone # Ascl2
# Trophoblast giant cells ::: features = c("Hand1","Stra13","Mdfi") , min.cutoff = 0)
# Labyrinthine trophoblast ::: features = c("Gcm1","Cebpa","Tfeb","Junb","Fosl1","Arnt", "Hif1a", "Hif1b", "Dlx3", "Pparg", "Pparbp") )


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


#write.csv(regulonTargetsInfo, paste0( Project, "_", "regulonTargetsInfo_ALL.csv"))
#write.csv(tableSubset, paste0( dataset, "_", "regulonTargetsInfo_interestingTFs.csv"))



message("--------------------------------------------------------------------------------")
message("+                           show TF expression                                  ")
message("+-------------------------------------------------------------------------------")

tsne_umap <- readRDS("int/tSNE_AUC_50pcs_50perpl.Rds")
colnames(tsne_umap$Y) <- c("tsne1" ,"tsne2")
tsne_SCENIC <- as.matrix(tsne_umap$Y)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
as.data.frame(rownames(aucell_regulonAUC))


# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tsne_SCENIC, exprMat, aucell_regulonAUC[rownames(aucell_regulonAUC) %in% c("Hif1a_extended (23g)", "Elf5_extended (1553g)", "Eomes (421g)", "Hand1 (12g)",  "Ets2_extended (124g)", "Junb (119g)"),], plots="Expression")


AUCell::AUCell_plotTSNE(tsne_SCENIC, exprMat, aucell_regulonAUC[rownames(aucell_regulonAUC) %in% c("Taf1 (632g)", "Rbpj (63g)", "Sox2_extended (19g)", "Kdm5a (185g)","Stat1 (157g)", "Irf7 (156g)"),], plots="Expression")

AUCell::AUCell_plotTSNE(tsne_SCENIC, exprMat, aucell_regulonAUC[rownames(aucell_regulonAUC) %in% c("Batf3_extended (90g)", "Sall4 (11g)", "Vezf1_extended (15g)", "Tcf7l1_extended (10g)","Gata3_extended (21g)", "Bhlhe40 (87g)"),], plots="Expression")




# plot AUC":::
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tsne_SCENIC, exprMat, aucell_regulonAUC[rownames(aucell_regulonAUC) %in% c("Hif1a_extended (23g)", "Elf5_extended (1553g)", "Eomes (421g)", "Hand1 (12g)",  "Ets2_extended (124g)", "Junb (119g)"),], plots="AUC")

AUCell::AUCell_plotTSNE(tsne_SCENIC, exprMat, aucell_regulonAUC[rownames(aucell_regulonAUC) %in% c("Taf1 (632g)", "Rbpj (63g)", "Sox2_extended (19g)", "Kdm5a (185g)","Stat1 (157g)", "Irf7 (156g)"),], plots="AUC")


AUCell::AUCell_plotTSNE(tsne_SCENIC, exprMat, aucell_regulonAUC[rownames(aucell_regulonAUC) %in% c("Batf3_extended (90g)", "Sall4 (11g)", "Vezf1_extended (15g)", "Tcf7l1_extended (10g)","Gata3_extended (21g)", "Bhlhe40 (87g)"),], plots="AUC")







message("--------------------------------------------------------------------------------")
message("+                     cluster specific regulons                                 ")
message("+-------------------------------------------------------------------------------")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- SCENIC::calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
rssPlot <- plotRSS(rss)
## Showing regulons and cell types with any RSS > 0.01 (dim: 6x5)
plotly::ggplotly(rssPlot$plot)



message("--------------------------------------------------------------------------------")
message("+                     investigate and Plot regulons                             ")
message("+-------------------------------------------------------------------------------")

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
regulon_names <- as.data.frame(cbind(onlyNonDuplicatedExtended(names(regulons))))
#regulon_names2 <- as.data.frame(cbind(names(regulons)))


#regulons <- loadInt(scenicOptions, "regulons")
names(regulons)
regulons[c("Elf5_extended", "Elf3","Gata1", "Ets2", "Junb")] # "Gcm1","Esrrb"

regulonNames <- c( "Stat2", "Cdx2", "Fos")

pdf(paste( "UMAP_regulonAUC_SCENIC__ORIG_together", "red.Eomes__green.Baf3__blue.Mta3",".pdf", sep="_"), width=15,height=11) # type= "cairo"
par(bg=NA)
regulonNames <- c( "Eomes", "Batf3","Mta3")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#text(-20, 5, attr(cellCol,"red"), col="red", cex=1.2, pos=4)
#text(8,-10, attr(cellCol,"green"), col="green3", cex=1.2, pos=4)
#text(3,8, attr(cellCol,"blue"), col="blue", cex=1.2, pos=4)
dev.off()


pdf(paste( "UMAP_regulonAUC_SCENIC__ORIG_together", "red.Eomes__green.Gata3__blue.Bhlhe40",".pdf", sep="_"), width=15,height=11)
par(bg=NA)
regulonNames <- c( "Eomes", "Gata3", "Bhlhe40") # Pou3f1 Elf5
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#text(-20, 5, attr(cellCol,"red"), col="red", cex=1.2, pos=4)
#text(-2,2, attr(cellCol,"blue"), col="blue", cex=1.2, pos=4)
#text(4,0, attr(cellCol,"green"), col="green3", cex=1.2, pos=4)
dev.off()


pdf(paste( "UMAP_regulonAUC_SCENIC__ORIG_together", "red_EomesSox2Zfp553Elf5Puo3f1__green_Elf2Stdb1Myc__blue_Foxo4Hif1aMsx2Mef2aSrebf2Mta3Rest",".pdf", sep="_"), width=15,height=11)
par(bg=NA)
regulonNames <- list(red=c("Eomes","Sox2", "Zfp553","Elf5","Pou3f1"),
                     green=c("Elf2", "Setdb1", "Myc"),
                     blue=c("Foxo4", "", "Hif1a", "Msx2", "Mef2a","Srebf2","Mta3","Rest"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
#text(-20, 5, attr(cellCol,"red"), col="red", cex=1.2, pos=4)
#text(-2,2, attr(cellCol,"blue"), col="blue", cex=1.2, pos=4)
#text(4,0, attr(cellCol,"green"), col="green3", cex=1.2, pos=4)
dev.off()







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



regulonNames <- list(red=c("Eomes", "Cdx2", "Sox2", "Tead4", "Elf5"),
                     green=c("Ascl2", "Ets2", "Gata1", "Batf3"),
                     blue=c("Gcm1", "Ovol2","Junb", "Hif1a", "Pparg", "Msx2"))

png(paste( "UMAP_regulonAUC.BINARY", dataset, "red.genes_TSC_green.genes.Cdx2EomesTead4Sox2_EPC.Ascl2__blue.genes_LP.Gcm1Ovol2", "v4",".png", sep="_"), width=2000,height=1500, type= "cairo")
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
text(-8, 0, attr(cellCol,"red"), col="black", cex=2.5, pos=4)
text(1,2, attr(cellCol,"green"), col="black", cex=2.5, pos=4)
text(4,3.5, attr(cellCol,"blue"), col="black", cex=2.5, pos=4)
dev.off()


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



regulonNames <- list(red=c("Cdx2","Tead4","Sox2"),
                     green=c("Ascl2"),
                     blue=c("Gcm1", "Ovol2"))
png(paste( "UMAP_regulonAUC", dataset, "red.genes_TSC_green.genes.Cdx2EomesTead4Sox2_EPC.Ascl2__blue.genes_LP.Gcm1Ovol2", "v4",".png", sep="_"), width=2000,height=1500, type= "cairo")
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.5)
text(-8, 0, attr(cellCol,"red"), col="black", cex=2.5, pos=4)
text(1,2, attr(cellCol,"green"), col="black", cex=2.5, pos=4)
text(4,3.5, attr(cellCol,"blue"), col="black", cex=2.5, pos=4)
dev.off()


# R only :::
regulonNames <- list(red=c("Cdx2","Eomes","Sox2"),
                     green=c("Vezf1","Gata1","Klf7"),
                     blue=c("Pparg", "Ghrl1")) # blue "Usf2","Nfya","Yy1","Maff"
png(paste( "UMAP_regulonAUC", dataset, "red.genes_TSC_green.genes.Cdx2EomesTead4Sox2_EPC.Ascl2__blue.genes_LP.Gcm1Ovol2", "v4",".png", sep="_"), width=2000,height=1500, type= "cairo")
par(bg=NA)
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.5)
text(-8, 0, attr(cellCol,"red"), col="black", cex=2.5, pos=4)
text(1,2, attr(cellCol,"green"), col="black", cex=2.5, pos=4)
text(4,3.5, attr(cellCol,"blue"), col="black", cex=2.5, pos=4)
dev.off()





regulonNames <- list(red=c("Eomes", "Cdx2", "Sox2"),
                     green=c("Myc", "Elf2", "Smarca4", "Setdb1"),
                     blue=c("Msx2", "Mef2a","Foxo4", "Hif1a", "", ""))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.5)
text(-20, 0, attr(cellCol,"red"), col="red", cex=1.5, pos=4)
text(1,-11, attr(cellCol,"green"), col="green", cex=1.5, pos=4)
text(3,4, attr(cellCol,"blue"), col="blue", cex=1.5, pos=4)





library(Seurat)
dr_coords <- Embeddings(seuratObject, reduction="tsne")

tfs <- c("Sox10","Irf1","Sox9", "Dlx5")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")










### HEATMAP::: Regulators for clusters or known cell types

library(seriation)


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUCx <- regulonAUC
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$SCT_snn_res.0.8),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

head(regulonActivity_byCellType_Scaled)
rownames(regulonActivity_byCellType_Scaled)

genes_for_heatmap <- c("E2f8","Foxo4","Rest","Myc","Srebf2","Elf2",  "Stat2", "Msx2", "Mef2a", "Hif1a", "Mta3","Klf6", "Setdb1", "Jund") #"Atf4",

#genes_for_heatmap <- c("Foxo3","Sp3","E2f8","Tead3","Grhl1","Ppard","Foxo4","Hdac2","Smarca4","Nfkb2","Creb3l2","Pou2f1","Myc","Max","Srebf2","E2f7","Fos","Elf2","Bhlhe40", "Elf5", "Ddit3", "Irf7", "Gata1", "Atf4", "Stat2", "Msx2", "Mef2a", "Hif1a", "Mta3","Klf6", "Setdb1", "Jund","Pou3f1")

#genes_for_heatmap <- c("Gcm1 (17g)", "Cdx2_extended (37g)", "Elf5 (31g)", "Eomes (14g)", "Ascl2_extended (14g)", "Sox2 (18g)", "Tead4 (454g)", "Ovol2_extended (17g)", "Smarca4 (158g)", "E2f7 (911g)", "Vezf1_extended (14g)","Irf1 (160g)","Smad1_extended (12g)","Gata1 (21g)","Zfp64 (12g)","E2f8 (17g)","Bdp1 (21g)","Jund (33g)","Mybl2 (12g)","Yy1 (733g)", "Rora_extended (12g)", "Cux1_extended (45g)","Foxo3 (24g)","Nr2f6_extended (10g)","Fbxl19 (14g)","Arnt (40g)","Klf7_extended (30g)","Cnot3_extended (71g)","Deaf1 (37g)","Atf5_extended (70g)","Junb (119g)")
#genes_for_heatmap <- c("Gcm1 (17g)", "Cdx2_extended (37g)", "Cdx2 (13g)", "Elf5 (31g)", "Eomes (14g)", "Ascl2_extended (14g)", "Sox2 (18g)", "Sox2 (11g)","Tead4 (454g)", "Tead4 (21g)","Junb (119g)","Ovol2_extended (17g)","Bhlhe40_extended (2092g)","Bhlhe40 (111g)","Tcf12_extended (34g)", "Jun (140g)", "Ets2 (1131g)", "Pparg (14g)" ,"Fos (40g)","Junb (127g)", "Elf5 (1101g)" , "Ybx1_extended (497g)",  "Tead3_extended (17g)") #, "Nfya (158g)","Xbp1 (73g)") # , 

all_TFs <- as.data.frame(rownames(regulonActivity_byCellType_Scaled))
colnames(all_TFs)[1] <- "TF"
all_TFs$TF_short <- all_TFs$TF
all_TFs$TF_short <- gsub( " \\(.*", "", all_TFs$TF_short)
all_TFs$TF_short <- gsub( "_extended", "", all_TFs$TF_short)

genes_for_heatmap <-  all_TFs[all_TFs$TF_short %in% genes_for_heatmap,]$TF 
#genes_for_heatmap <-  all_TFs$TF

regulonActivity_byCellType_Scaled2 <- regulonActivity_byCellType_Scaled[rownames(regulonActivity_byCellType_Scaled) %in% genes_for_heatmap,]


o = seriate(max(regulonActivity_byCellType_Scaled2) - regulonActivity_byCellType_Scaled2, method = "BEA_TSP")
#Heatmap(max(mat) - mat, name = "mat", 
#        row_order = get_order(o, 1), column_order = get_order(o, 2),
#        column_title = "seriation by BEA_TSP method")


ht1 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled2, color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100), treeheight_row=10, treeheight_col=10, border_color=NA, row_order = get_order(o, 1), column_order = get_order(o, 2))

ht1

pdf(paste( "Pheatmap_regulonAUC", Project, "64k_subsData_","res.0.8_ALL_TFs",".pdf", sep="_"), width=10,height=18)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht1
dev.off()





# cluster 20 removed:::
regulonActivity_byCellType_Scaled3 <- regulonActivity_byCellType_Scaled2[,-10]

o = seriate(max(regulonActivity_byCellType_Scaled3) - regulonActivity_byCellType_Scaled3, method = "BEA_TSP")

ht2 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled3,  color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),  treeheight_row=10, treeheight_col=10, border_color=NA, column_order = get_order(o, 2), row_order = get_order(o, 1)) # 

ht2


pdf(paste( "Pheatmap_regulonAUC", Project, "64k_subsData_","res.0.8_ALL_TFs_clust20rm",".pdf", sep="_"), width=10,height=18)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht2
dev.off()


regulonActivity_byCellType_Scaled4 <- regulonActivity_byCellType_Scaled3[,c("2","7", "8","14","9","10","1","4","5","6","16","19","13")]
o = seriate(max(regulonActivity_byCellType_Scaled4) - regulonActivity_byCellType_Scaled4, method = "BEA_TSP")

ht2 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled4,  color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),  treeheight_row=10, treeheight_col=10, border_color=NA, column_order = get_order(o, 2), row_order = get_order(o, 1), angle_col = 0)  # 



pdf(paste( "Pheatmap_regulonAUC", Project, "64k_subsData_","res.0.8_ALL_TFs_clust20rm_EPC_vs_LP",".pdf", sep="_"), width=6,height=4.5)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht2
dev.off()





message("--------------------------------------------------------------------------------")
message("+                        plot binarised heatmap                                 ")
message("+-------------------------------------------------------------------------------")



topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators, "topRegulators_I_celltype.pred.csv")

#Binarized version (~ percentage of cells of that cell type/cluster with the regulon active)
minPerc <- 0.7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype.pred), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

ht2 <- pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5, 
                          color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                          treeheight_row=10, treeheight_col=10, border_color=NA)

dim(binaryActPerc_subset)

pdf(paste( "Pheatmap_regulonAUC_Binarized", "minPerc", minPerc, "celltype.pred", dataset, "Clusters1.25_",".pdf", sep="_"), width=6,height=8)
par(bg=NA)
ht2
dev.off()



#Binarized version (~ percentage of cells of that cell type/cluster with the regulon active)
minPerc <- 0.3
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$ClustIdent), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

ht2 <- pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5, 
                          color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                          treeheight_row=10, treeheight_col=10, border_color=NA)

dim(binaryActPerc_subset)

pdf(paste( "Pheatmap_regulonAUC_Binarized", "minPerc", minPerc, "ClustIdent", dataset, "Clusters1.25_",".pdf", sep="_"), width=12,height=15)
par(bg=NA)
ht2
dev.off()




message("--------------------------------------------------------------------------------")
message("+         Complex heatmap for genes within selected regulons                    ")
message("+-------------------------------------------------------------------------------")
#mat_clusterMeans <- AverageExpression(matrix.su, assays = "SCT", slot = "scale.data")
#saveRDS(mat_clusterMeans, "mat_clusterMeans_harmony.orig_500g_res0.8.Rds")

mat_clusterMeans <- readRDS("mat_clusterMeans_harmony.orig_500g_res0.8.Rds")

regulons_incidMat <- readRDS("int/2.6_regulons_asIncidMat.Rds")
regulons <- regulonsToGeneLists(regulons_incidMat)


regulonsSelected <- selectRegulons(regulons, c("Ets2", "Junb"))

mat_clusterMeans_SCT <- mat_clusterMeans[[1]][rownames(mat_clusterMeans[[1]]) %in% unique(unlist(regulonsSelected)),]
mat_clusterMeans_SCT[1:5, 1:5]
mat_clusterMeans_SCT <- mat_clusterMeans_SCT - rowMeans(mat_clusterMeans_SCT)
min_val <- min(mat_clusterMeans_SCT)
max_val <- max(mat_clusterMeans_SCT)

clust_order <-c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8") # for R 500g res0.8
#clust_order <-c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0") # for I 500g res0.8
mat_clusterMeans_SCT <- mat_clusterMeans_SCT[, match(clust_order, colnames(mat_clusterMeans_SCT))]

# for slot: data:
#f1 = colorRamp2( c(-10, -3, 0, 3, 10), c("#173715", "green3", "grey95", "#bf7acd", "#4b0e81"), space = "LAB") 
# for scale.data:
f1 = colorRamp2( c(min_val, 0, 3, 10), c("green3", "grey95", "#bf7acd", "#4b0e81"), space = "LAB") 
lgd1 = Legend(col_fun = f1, title = "Mean cluster expression", at = c(min_val, max_val )  )

ht1 = Heatmap((mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "Marker genes (abs l2fc > 0.6)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE ,  row_names_side ="left",  row_title_rot = 90, column_title_rot = 0 , na_col = "lightgrey", row_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "bottom", column_names_side = "top", row_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1
ht1

pdf(paste(Project, dataset_name, "Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_clusterOrderedTempora.pdf", sep="_"), onefile=FALSE, width=7, height=45) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


ht1_t = Heatmap(t(mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "Marker genes (abs l2fc > 0.6)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = FALSE ,  row_names_side ="left",  row_title_rot = 90, column_title_rot = 0 , na_col = "lightgrey", column_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "top", column_names_side = "bottom",  row_dend_reorder = TRUE, column_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1
ht1_t

pdf(paste(Project, dataset_name,"t_Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_clusterOrderedTempora2.pdf", sep="_"), onefile=FALSE, width=50, height=7) 
par(bg=NA)
draw(ht1_t, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()





message("--------------------------------------------------------------------------------")
message("+         Complex heatmap comparing active regulons between I and R             ")
message("+-------------------------------------------------------------------------------")

# steps... 
# 1. load matrix.su_sep for R or I
# 2. get regulon activity per cell - and subset by R or I
# 3. plot heatmap -R or I

library("ComplexHeatmap")
library("circlize")





message("--------------------------------------------------------------------------------")
message("+        AUCell activity by treatment...... HEATMAP            ")
message("+-------------------------------------------------------------------------------")

cellInfo$Clust_Treatment <- paste(cellInfo$SCT_snn_res.0.8, cellInfo$Treatment, sep = "_")
table(cellInfo$Clust_Treatment)

library(seriation)


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUCx <- regulonAUC
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Clust_Treatment),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

head(regulonActivity_byCellType_Scaled)
rownames(regulonActivity_byCellType_Scaled)

genes_for_heatmap <- c("E2f8","Foxo4","Rest","Myc","Srebf2","Elf2",  "Stat2", "Msx2", "Mef2a", "Hif1a", "Mta3","Klf6", "Setdb1", "Jund") #"Atf4",

#genes_for_heatmap <- c("Foxo3","Sp3","E2f8","Tead3","Grhl1","Ppard","Foxo4","Hdac2","Smarca4","Nfkb2","Creb3l2","Pou2f1","Myc","Max","Srebf2","E2f7","Fos","Elf2","Bhlhe40", "Elf5", "Ddit3", "Irf7", "Gata1", "Atf4", "Stat2", "Msx2", "Mef2a", "Hif1a", "Mta3","Klf6", "Setdb1", "Jund","Pou3f1")

#genes_for_heatmap <- c("Gcm1 (17g)", "Cdx2_extended (37g)", "Elf5 (31g)", "Eomes (14g)", "Ascl2_extended (14g)", "Sox2 (18g)", "Tead4 (454g)", "Ovol2_extended (17g)", "Smarca4 (158g)", "E2f7 (911g)", "Vezf1_extended (14g)","Irf1 (160g)","Smad1_extended (12g)","Gata1 (21g)","Zfp64 (12g)","E2f8 (17g)","Bdp1 (21g)","Jund (33g)","Mybl2 (12g)","Yy1 (733g)", "Rora_extended (12g)", "Cux1_extended (45g)","Foxo3 (24g)","Nr2f6_extended (10g)","Fbxl19 (14g)","Arnt (40g)","Klf7_extended (30g)","Cnot3_extended (71g)","Deaf1 (37g)","Atf5_extended (70g)","Junb (119g)")
#genes_for_heatmap <- c("Gcm1 (17g)", "Cdx2_extended (37g)", "Cdx2 (13g)", "Elf5 (31g)", "Eomes (14g)", "Ascl2_extended (14g)", "Sox2 (18g)", "Sox2 (11g)","Tead4 (454g)", "Tead4 (21g)","Junb (119g)","Ovol2_extended (17g)","Bhlhe40_extended (2092g)","Bhlhe40 (111g)","Tcf12_extended (34g)", "Jun (140g)", "Ets2 (1131g)", "Pparg (14g)" ,"Fos (40g)","Junb (127g)", "Elf5 (1101g)" , "Ybx1_extended (497g)",  "Tead3_extended (17g)") #, "Nfya (158g)","Xbp1 (73g)") # , 

all_TFs <- as.data.frame(rownames(regulonActivity_byCellType_Scaled))
colnames(all_TFs)[1] <- "TF"
all_TFs$TF_short <- all_TFs$TF
all_TFs$TF_short <- gsub( " \\(.*", "", all_TFs$TF_short)
all_TFs$TF_short <- gsub( "_extended", "", all_TFs$TF_short)

genes_for_heatmap <-  all_TFs[all_TFs$TF_short %in% genes_for_heatmap,]$TF 
#genes_for_heatmap <-  all_TFs$TF

regulonActivity_byCellType_Scaled2 <- regulonActivity_byCellType_Scaled[rownames(regulonActivity_byCellType_Scaled) %in% genes_for_heatmap,]


o = seriate(max(regulonActivity_byCellType_Scaled2) - regulonActivity_byCellType_Scaled2, method = "BEA_TSP")
#Heatmap(max(mat) - mat, name = "mat", 
#        row_order = get_order(o, 1), column_order = get_order(o, 2),
#        column_title = "seriation by BEA_TSP method")


ht1 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled2, color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100), treeheight_row=10, treeheight_col=10, border_color=NA, row_order = get_order(o, 1), column_order = get_order(o, 2))

ht1

pdf(paste( "Pheatmap_regulonAUC", Project, "64k_subsData_","res.0.8_ALL_TFs",".pdf", sep="_"), width=10,height=18)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht1
dev.off()





# cluster 20 removed:::
regulonActivity_byCellType_Scaled3 <- regulonActivity_byCellType_Scaled2[,-10]

o = seriate(max(regulonActivity_byCellType_Scaled3) - regulonActivity_byCellType_Scaled3, method = "BEA_TSP")

ht2 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled3,  color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),  treeheight_row=10, treeheight_col=10, border_color=NA, column_order = get_order(o, 2), row_order = get_order(o, 1)) # 

ht2


pdf(paste( "Pheatmap_regulonAUC", Project, "64k_subsData_","res.0.8_ALL_TFs_clust20rm",".pdf", sep="_"), width=10,height=18)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht2
dev.off()


regulonActivity_byCellType_Scaled4 <- regulonActivity_byCellType_Scaled3[,c("2","7", "8","14","9","10","1","4","5","6","16","19","13")]
o = seriate(max(regulonActivity_byCellType_Scaled4) - regulonActivity_byCellType_Scaled4, method = "BEA_TSP")

ht2 <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled4,  color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),  treeheight_row=10, treeheight_col=10, border_color=NA, column_order = get_order(o, 2), row_order = get_order(o, 1), angle_col = 0)  # 



pdf(paste( "Pheatmap_regulonAUC", Project, "64k_subsData_","res.0.8_ALL_TFs_clust20rm_EPC_vs_LP",".pdf", sep="_"), width=6,height=4.5)
#png(paste( "Pheatmap_regulonAUC", dataset, "ClustIdent","Clusters1.25_interestingTFs",".png", sep="_"), width=900,height=700, type = "cairo")
par(bg=NA)
ht2
dev.off()






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
