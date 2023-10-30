#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())

############  this script makes matrix from realigned batches 
############  (read length unified to 50bp to remove batch effect): 
############  Dropseq 0004, 0006, 0007, 0008, 0011, 0014. 


message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
  library("dplyr")
  library("methods")
  library("utils")
  library("ggplot2")
  library("ggrepel")
  library("cowplot")
  library("Matrix")
  library("matrixStats")
  library("Seurat")
  library("parallel")
  library("devtools")
  library("grDevices")
  library("ggalt")
  library("Matrix.utils")
  #library("biomaRt")
})


genecut       <- 300
mincells      <- 3
resolution    <- "1"
normalisation <- "log2"

cores <- detectCores()
print(cores)


#baseDir <- "/home/mn367/rds/rds-dscj1-dropseq"
baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
print(baseDir)
setwd(baseDir)



message("--------------------------------------------------------------------------------")
message("+                                Load objects:::                                ")
message("+-------------------------------------------------------------------------------")

matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_40pc_alignRate_53.Rds")
sampleTable <- read.csv("sampleTable_MERGED_53.csv")


message("--------------------------------------------------------------------------------")
message("+                          Sample files on ctr-bfx                              ")
message("+-------------------------------------------------------------------------------")

sampleTable <- read.csv("sampleTable_183samples_all_batches_alignrate.csv")
sampleTable <- sampleTable[,-10]

sampleTable$Batch2 <- ifelse(sampleTable$Batch == "Batch0006" | sampleTable$Batch == "Batch0008", "A00", "000")
sampleTable$Batch2 <- ifelse(sampleTable$Batch == "Batch0007" | sampleTable$Batch == "Batch0011", "0B0", sampleTable$Batch2)
sampleTable$Batch2 <- ifelse(sampleTable$Batch == "Batch0004" | sampleTable$Batch == "Batch0014", "00C", sampleTable$Batch2)

merged_samples_df <- read.csv("Final_set_to_merge_MERGED.csv")
merged_samples_molten <- reshape2::melt(merged_samples_df, id.vars = c("Name_of_rep_set", "Merged.File"))
merged_samples_molten <- na.omit(merged_samples_molten)

sampleTable$MergedRepSet <- merged_samples_molten[match(sampleTable$sampleLabels, merged_samples_molten$value),]$Name_of_rep_set
sampleTable$MergedFilePath <- merged_samples_molten[match(sampleTable$sampleLabels, merged_samples_molten$value),]$Merged.File
sampleTable$MergedFilePath2 <- ifelse(is.na(sampleTable$MergedFilePath), as.character(sampleTable$fileNameDGE), as.character(sampleTable$MergedFilePath))
sampleTable$Nxxx <- gsub(".*_N", "N" , sampleTable$sampleLabels)
sampleTable$sampleLabels2 <-  ifelse(is.na(sampleTable$MergedRepSet), as.character(sampleTable$sampleLabels), paste(as.character(sampleTable$MergedRepSet),sampleTable$Nxxx, sep = "_"))

change_batch <- sampleTable[,c("MergedRepSet","Batch", "Batch2")]
change_batch$Batch3 <- change_batch$Batch2
change_batch$Batch3[change_batch$MergedRepSet== "Dups_1Aa_reps"] <- "A0C"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_1Ab_reps"] <- "A00"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_setB_reps"] <- "A00"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_setA1_reps"] <- "0BC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_1C_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_1D_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_setA2_reps"] <- "A00"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_1B_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_2C_reps"] <- "AB0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_2D_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_2B_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_5A_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_5C_reps"] <- "AB0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_6C_reps"] <- "AB0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_7A_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_6B_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_8A_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_7B_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_8D_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_6A_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_5B_reps"] <- "A0C"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_6D_reps"] <- "A0C"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_8C_reps"] <- "AB0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_8B_reps"] <- "ABC"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48R3_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48I2_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48R2_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36R2_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36I4_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36I3_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48R1_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48I4_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48I3_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36R4_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36R1_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36I1_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36R3_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48R4_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_36I2_reps"] <- "0B0"
change_batch$Batch3[change_batch$MergedRepSet== "Dups_48I1_reps"] <- "0B0"
change_batch <- change_batch[!duplicated(change_batch$MergedRepSet),]
change_batch <- na.omit(change_batch)                   

sampleTable_MERGED <- sampleTable[!duplicated(sampleTable$sampleLabels2),]
sampleTable_MERGED$Merged_batch <- change_batch[match(sampleTable_MERGED$MergedRepSet, change_batch$MergedRepSet),]$Batch3
sampleTable_MERGED$Merged_batch <- ifelse(is.na(sampleTable_MERGED$Merged_batch), as.character(sampleTable_MERGED$Batch2) , sampleTable_MERGED$Merged_batch)
sampleTable_MERGED <- sampleTable_MERGED[,c("sampleLabels2","CellType","Age","Treatment","MergedRepSet", "MergedFilePath2", "Merged_batch" , "Nxxx"  )]


sampleTable_MERGED$MergedFilePath <- gsub( ".*\\.N700\\/", "\\/storage\\/CTR-Projects\\/CTR_DropSeq\\/NewMegatron\\/NewAlignments_50bp\\/ForMatrix\\/Merged_samples\\/", sampleTable_MERGED$MergedFilePath2)



message("--------------------------------------------------------------------------------")
message("+                          Make matrix now                                      ")
message("+-------------------------------------------------------------------------------")


head(sampleTable,3)

sampleFiles <- sampleTable$MergedFilePath
length(sampleFiles)
names(sampleFiles) <- sampleTable$sampleLabels2
sampleNames <- sampleTable$sampleLabels2



Sample_Cell_Genes.tbl <- data.frame(sampleName=character(),Genes=integer(),Cells=integer() )

matrix.su <- ""
count     <- 0


#file <- sampleFiles[[2]]
#sampleName <- sampleNames[[2]]

for(file in sampleFiles)
{
  sampleName <- names(sampleFiles)[sampleFiles == file]
  
  message(paste0("-- Running ", count, " : ", sampleName))
  tmp.tab <- read.table(gzfile(paste0(file)), sep="\t",header=TRUE, row.names=1, stringsAsFactors = TRUE)
  
  message( paste0("-- Genes=", dim(tmp.tab)[1], " Cells=", dim(tmp.tab)[2]),
           " FilteredCells=", sum(colSums(tmp.tab) > genecut))
  
  tmp.cellgene                    <- ""
  tmp.cellgene                    <- data.frame(sampleName, dim(tmp.tab)[1], dim(tmp.tab)[2])
  colnames(tmp.cellgene)          <- c("sampleName", "Genes", "Cells")
  Sample_Cell_Genes.tbl           <- rbind(Sample_Cell_Genes.tbl, tmp.cellgene)
  colnames(Sample_Cell_Genes.tbl) <- c("sampleName", "Genes", "Cells")
  
  if(count == 0)
  {
    message("-- Creating matrix.su")
    matrix.su <- ""
    matrix.su <- CreateSeuratObject(counts = tmp.tab, project = sampleName, min.cells = mincells, min.features = genecut) # not `counts` but `raw.data` in seurat2.3.4
    matrix.su <- RenameCells(matrix.su, add.cell.id = sampleName)
    message( paste0("-- SU Size= ", format( object.size(matrix.su), units='auto') ))
  }
  else
  {
    message("-- Appending matrix.su")
    
    merged.su <- CreateSeuratObject(counts = tmp.tab, min.cells = mincells, min.features = genecut, project = sampleName)
    merged.su <- RenameCells(merged.su, add.cell.id = sampleName)
    matrix.su <- merge(x=matrix.su, y=merged.su, project = "Troph")
    #matrix.su <- MergeSeurat(matrix.su, merged.su, project = "Troph", add.cell.id2 = c(sampleName))
    message( paste0("-- SU Size= ", format( object.size(matrix.su), units='auto') ))
  }
  
  tmp.tab <- ""
  count = count+1
}



message("-- Seurat matrix done")




print(count)
print( table(matrix.su@meta.data$orig.ident) )


message("--------------------------------------------------------------------------------")
message("+       Annotate  matrix.su                               ")
message("+-------------------------------------------------------------------------------")

matrix.su@meta.data$orig.ident2 <- gsub( "Dropseq", "Batch", matrix.su@meta.data$orig.ident)
matrix.su@meta.data$Experiment <- sampleTable[match(matrix.su@meta.data$orig.ident, sampleTable$sampleLabels2),]$CellType
matrix.su@meta.data$Age <- sampleTable[match(matrix.su@meta.data$orig.ident, sampleTable$sampleLabels2),]$Age
matrix.su@meta.data$Treatment <- sampleTable[match(matrix.su@meta.data$orig.ident, sampleTable$sampleLabels2),]$Treatment
matrix.su@meta.data$Batch <- sampleTable[match(matrix.su@meta.data$orig.ident, sampleTable$sampleLabels2),]$Merged_batch


head(matrix.su@meta.data,2)
unique(matrix.su@meta.data$Batch)
table(matrix.su@meta.data[,c(5,8)])


#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_40pc_alignRate_53.Rds")
#write.csv(Sample_Cell_Genes.tbl, "Sample_Cell_Genes.tbl__MergedSamples_40pc_alignRate_53.csv")



message("--------------------------------------------------------------------------------")
message("+       SUBSETTING                               ")
message("+-------------------------------------------------------------------------------")

nrow(matrix.su@meta.data) #  230037
min(matrix.su@meta.data$nCount_RNA) # 294
min(matrix.su@meta.data$nFeature_RNA) #  272

keep_cells <- rownames(matrix.su@meta.data[ matrix.su@meta.data$nFeature_RNA >= 300 ,])
length(keep_cells)  #  229793
matrix.su <- SubsetData(matrix.su, cells = keep_cells)
nrow(matrix.su@meta.data) #  

#saveRDS(matrix.su, "CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53.Rds")



