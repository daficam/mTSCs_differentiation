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

matrix.su <- readRDS("CTR_dscj1_Realigned_50bp_matrix.su.Rds")

sampleTable <- read.csv("sampleTable2.csv")



message("--------------------------------------------------------------------------------")
message("+                          Sample files on ctr-bfx                              ")
message("+-------------------------------------------------------------------------------")

sampleTable <- read.csv("sampleTable.csv")
sampleTable$fileNameDGE_ctrbfx <- gsub("\\/home\\/mn367\\/rds\\/rds-dscj1-dropseq\\/CTR_DropSeq", "\\/storage\\/CTR-Projects\\/CTR_DropSeq\\/NewMegatron\\/NewAlignments_50bp\\/ForMatrix\\/Dropseq", sampleTable$fileNameDGE)
sampleTable$fileNameDGE_ctrbfx <- gsub("SLX.*", "", sampleTable$fileNameDGE_ctrbfx)
sampleTable$fileNameDGE_ctrbfx <- paste0(sampleTable$fileNameDGE_ctrbfx, sampleTable$fileName_no_path)

#sampleTable$fileName_no_path <- gsub( ".*\\/", "", sampleTable$fileNameDGE)

#sampleFiles_df <- as.data.frame(list.files(baseDir, pattern='*N7[0-9][0-9].final.dge.txt.gz', recursive = TRUE))
#colnames(sampleFiles_df) <- "sampleFiles_ctrbfx"
#sampleFiles_df$sampleFiles_no_path <- gsub( ".*\\/", "", sampleFiles_df$sampleFiles_ctrbfx)

#nrow(sampleFiles_df[sampleFiles_df$sampleFiles_no_path %in% sampleTable$fileName_no_path,])
#sampleFiles_df <- sampleFiles_df[sampleFiles_df$sampleFiles_no_path %in% sampleTable$fileName_no_path,]

#sampleTable$fileNameDGE_ctrbfx <- sampleFiles_df[match(sampleTable$fileName_no_path,sampleFiles_df$sampleFiles_no_path ),]$sampleFiles_ctrbfx

sampleFiles <- sampleTable$fileNameDGE_ctrbfx
length(sampleFiles)


message("--------------------------------------------------------------------------------")
message("+                          Sample files on ctr-bfx                              ")
message("+-------------------------------------------------------------------------------")

samples_Dropseq_0008       <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/ForMatrix/Dropseq_0008"
sampleFiles_Dropseq_0008 <- list.files(samples_Dropseq_0008, pattern='*N7[0-9][0-9].final.dge.txt.gz', recursive = TRUE)
sampleFiles_Dropseq_0008 <- paste0(samples_Dropseq_0008, "/", sampleFiles_Dropseq_0008)

samples_Dropseq_0007       <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/ForMatrix/Dropseq_0007"
sampleFiles_Dropseq_0007 <- list.files(samples_Dropseq_0007, pattern='*N7[0-9][0-9].final.dge.txt.gz', recursive = TRUE)
sampleFiles_Dropseq_0007 <- paste0(samples_Dropseq_0007, "/", sampleFiles_Dropseq_0007)

samples_Dropseq_0006       <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/ForMatrix/Dropseq_0006"
sampleFiles_Dropseq_0006 <- list.files(samples_Dropseq_0006, pattern='*N7[0-9][0-9].final.dge.txt.gz', recursive = TRUE)
sampleFiles_Dropseq_0006 <- paste0(samples_Dropseq_0006, "/", sampleFiles_Dropseq_0006)

samples_Dropseq_0004       <-  "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/ForMatrix/Dropseq_0004"
sampleFiles_Dropseq_0004 <- list.files(samples_Dropseq_0004, pattern='*N7[0-9][0-9].final.dge.txt.gz', recursive = TRUE)
sampleFiles_Dropseq_0004 <- paste0(samples_Dropseq_0004, "/", sampleFiles_Dropseq_0004)

samples_Dropseq_0014      <-  "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/ForMatrix/Dropseq_0014"
sampleFiles_Dropseq_0014 <- list.files(samples_Dropseq_0014, pattern='*N7[0-9][0-9].final.dge.txt.gz', recursive = TRUE)
sampleFiles_Dropseq_0014 <- paste0(samples_Dropseq_0014, "/", sampleFiles_Dropseq_0014)

samples_Dropseq_0011      <-  "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/ForMatrix/Dropseq_0011"
sampleFiles_Dropseq_0011 <- list.files(samples_Dropseq_0011, pattern='*N7[0-9][0-9].final.dge.txt.gz', recursive = TRUE)
sampleFiles_Dropseq_0011 <- paste0(samples_Dropseq_0011, "/", sampleFiles_Dropseq_0011)


sampleFiles <- unique(c(sampleFiles_Dropseq_0011,
                        sampleFiles_Dropseq_0008, 
                        sampleFiles_Dropseq_0007,
                        sampleFiles_Dropseq_0006,
                        sampleFiles_Dropseq_0004,
                        sampleFiles_Dropseq_0014))
# 190 now!!!



sampleNames <- gsub(".final.dge.txt.gz", "", sampleFiles)
sampleNames <- gsub("\\/storage\\/CTR-Projects\\/CTR_DropSeq\\/NewMegatron\\/NewAlignments_50bp\\/ForMatrix\\/", "", sampleNames)
sampleNames <- gsub("\\/SLX-", "_", sampleNames)
sampleNames <- gsub("Dropseq_", "Batch", sampleNames)
sampleNames <- gsub(".[ATXGC]{8}.", "_", sampleNames)

head(sampleNames, 40)
length(sampleNames)
length(unique(sampleNames))


sampleTable <- data.frame(sampleLabels=sampleNames, fileNameDGE=sampleFiles)
head(sampleTable)

sampleTable$Batch <- sampleTable$sampleLabels
sampleTable$Batch <- gsub("_.*","",      sampleTable$Batch)

sampleTable$CellType <- sampleTable$sampleLabels
sampleTable$CellType <- gsub("Batch[0-9]{4}_","",      sampleTable$CellType)

sampleTable$CellType <- gsub("15980_N707","R.24",      sampleTable$CellType)
sampleTable$CellType <- gsub("15980_N722","I.24",      sampleTable$CellType)
sampleTable$CellType <- gsub("15980_N723","R.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15980_N724","R.1",       sampleTable$CellType)
sampleTable$CellType <- gsub("15980_N726","R.4",       sampleTable$CellType)
sampleTable$CellType <- gsub("15980_N727","15980_xN727", sampleTable$CellType)
sampleTable$CellType <- gsub("15980_N729","R.24",      sampleTable$CellType)
sampleTable$CellType <- gsub("15982_N712","R.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15982_N716","I.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15982_N720","R.4",       sampleTable$CellType)
sampleTable$CellType <- gsub("15982_N722","R.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15982_N723","I.1",       sampleTable$CellType)
sampleTable$CellType <- gsub("15982_N726","R.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15982_N728","I.24",      sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N701","I.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N703","I.1",       sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N704","I.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N705","I.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N706","R.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N707","R.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N729","R.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15985_N714","LNA.LNA",   sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N702","I.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N703","I.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N705","I.4",      sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N706","I.24",      sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N710","R.36",      sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N711","R.1",       sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N715","15987_xN715", sampleTable$CellType)
sampleTable$CellType <- gsub("15987_N718","LNA.LNA",   sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N704","0.0",       sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N710","I.24",      sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N712","I.1",       sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N716","R.1",       sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N718","I.1",       sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N719","I.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N720","I.48",      sampleTable$CellType)
sampleTable$CellType <- gsub("15988_N721","LNA.LNA",   sampleTable$CellType)


sampleTable$CellType <- gsub("1039_N700","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("1030_N715","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("1029.N727","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("1004.N716","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("1006.N712","I.1",        sampleTable$CellType)

# check if 9617 == 9340
# samples : N711, N703, N714 - ???
sampleTable$CellType <- gsub("9617_N718","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N719","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N720","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N721","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N722","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N724","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N726","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N727","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("9617_N729","0.0",        sampleTable$CellType)

sampleTable$CellType <- gsub("7635_N723","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N712","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N703","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N718","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N702","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N711","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N724","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N726","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N719","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N704","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N701","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N727","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N705","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N726","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N715","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N720","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N728","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N710","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N706","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N722","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N729","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N707","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N723","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N712","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N703","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N718","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N702","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N711","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N724","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N716","R.1",        sampleTable$CellType)

sampleTable$CellType <- gsub("7633_N723","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N712","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N703","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N718","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N702","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N711","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N724","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N716","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N719","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N704","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N701","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N705","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N726","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N715","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N720","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N728","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N710","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N706","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N722","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N729","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N707","R.24",        sampleTable$CellType)

sampleTable$CellType <- gsub("7637_N711","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N718","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N728","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N724","0.0",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N701","I.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N702","I.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N703","I.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N704","I.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N706","R.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N707","R.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N710","R.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N712","R.36",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N705","I.48",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N716","I.48",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N719","I.48",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N720","I.48",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N722","R.48",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N723","R.48",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N726","R.48",        sampleTable$CellType)
sampleTable$CellType <- gsub("7637_N729","R.48",        sampleTable$CellType)

sampleTable$CellType <- gsub("9339_N723","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N712","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N718","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N711","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N724","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N716","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N727","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N714","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N719","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N726","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N715","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N720","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N728","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N710","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N722","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N729","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N707","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N721","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9339_N701","0.0",        sampleTable$CellType)


sampleTable$CellType <- gsub("9616_N723","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N712","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N718","I.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N711","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N724","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N716","R.1",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N727","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N714","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N719","I.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N726","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N715","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N720","R.4",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N728","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N710","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N722","I.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N729","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N707","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N721","R.24",        sampleTable$CellType)
sampleTable$CellType <- gsub("9616_N701","0.0",        sampleTable$CellType)

#sampleTable$CellType <- gsub("x704.x704","0.0",        sampleTable$CellType)
#sampleTable$CellType <- gsub("x715.x715", NA,        sampleTable$CellType)
#sampleTable$CellType <- gsub("x727.x727", NA,        sampleTable$CellType)


sampleTable$CellType <- gsub("10298_N.*", NA,        sampleTable$CellType)
sampleTable$CellType <- gsub("7633_N.*", NA,        sampleTable$CellType)
sampleTable$CellType <- gsub("7635_N.*",NA ,        sampleTable$CellType) ##
sampleTable$CellType <- gsub("9617_N.*",NA,        sampleTable$CellType)


sampleTable$CellType <- gsub("15305_N701", "0.0",        sampleTable$CellType)


sampleTable <- sampleTable[!is.na(sampleTable$CellType),] # 181 samples now


Age.Treat             <- t(data.frame(strsplit(sampleTable$CellType, ".", fixed=TRUE)))
colnames(Age.Treat)   <- c("Treatment", "Age")
Age.Treat             <- data.frame(Age.Treat)
sampleTable$Age       <- Age.Treat$Age
sampleTable$Treatment <- Age.Treat$Treatment
#head(sampleTable, 10)
nrow(sampleTable)
rownames(sampleTable) <- sampleTable$sampleLabels
head(sampleTable)
#print(sampleTable)

sampleTable$sampleName <- gsub("Batch[0-9]{4}_","", sampleTable$sampleLabels)

#write.csv(sampleTable, "sampleTable_183samples_all_batches.csv")
#write.csv(sampleTable[sampleTable$Batch != "Batch0004",], "sampleTable_172_samples_0004rm.csv")

#sampleTable <- subset(sampleTable, sampleTable$Batch != "Batch0004")

#sampleTable <- subset(sampleTable, sampleTable$Batch == "Batch0011")
#sampleTable <- subset(sampleTable, sampleTable$Treatment == "R")


sampleTable <- read.csv("sampleTable.csv")
sampleFiles <- sampleTable$fileNameDGE
length(sampleFiles)
sampleNames <- sampleTable$sampleLabels
length(sampleNames)


#file <- sampleFiles[[1]]
message("-- Sample table done")

#tmp.tab <- read.table(gzfile(paste0(file)), sep="\t",header=TRUE, row.names=1, stringsAsFactors = TRUE)











message("--------------------------------------------------------------------------------")
message("+                          Make matrix now                                      ")
message("+-------------------------------------------------------------------------------")



Sample_Cell_Genes.tbl <- data.frame(sampleName=character(),Genes=integer(),Cells=integer() )

matrix.su <- ""
count     <- 0


names(sampleFiles) <- sampleTable$sampleName


file <- sampleFiles[[2]]

for(file in sampleFiles)
{
  sampleName <- file
  sampleName <- gsub(".final.dge.txt.gz", "", sampleName)
  #sampleName <- gsub("\\/home\\/mn367\\/rds\\/rds-dscj1-dropseq\\/", "", sampleName)
  sampleName <- gsub("\\/storage\\/CTR-Projects\\/CTR_DropSeq\\/NewMegatron\\/NewAlignments_50bp\\/ForMatrix\\/", "", sampleName)
  sampleName <- gsub("\\/SLX-", "_", sampleName)
  sampleName <- gsub("Dropseq_", "Batch", sampleName)
  sampleName <- gsub(".[ATXGC]{8}.", "_", sampleName)
  
  
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








setwd(resDir)


print(count)
print( table(matrix.su@meta.data$orig.ident) )

#¢saveRDS(matrix.su, file = paste(resDir, "/matrix.su__181_samples_seurat.2.3.4__400g.rds", sep = ""))
#saveRDS(matrix.su, file = paste(resDir, "/matrix.su__181_samples_seurat.2.3.4__300g.rds", sep = ""))
#slotNames(matrix.su)


saveRDS(matrix.su, "CTR_dscj1_Realigned_50bp_matrix.su_183samples.Rds")
write.csv(Sample_Cell_Genes.tbl, "Sample_Cell_Genes.tbl_181_samples_x.csv")

#rownames(matrix.su@meta.data) <- gsub( "\\/storage\\/CTR-Projects\\/CTR_DropSeq\\/NewMegatron\\/NewAlignments_50bp\\/", "", rownames(matrix.su@meta.data))
#matrix.su@meta.data$orig.ident <- gsub( "\\/storage\\/CTR-Projects\\/CTR_DropSeq\\/NewMegatron\\/NewAlignments_50bp\\/", "", matrix.su@meta.data$orig.ident)




message("--------------------------------------------------------------------------------")
message("+       SUBSETTING                               ")
message("+-------------------------------------------------------------------------------")


##### below is just subsetting, done now so do not run:::

matrix.su@meta.data$orig.ident2 <- gsub( "Dropseq", "Batch", matrix.su@meta.data$orig.ident)

matrix.su@meta.data$orig.ident2 <- gsub( "Dropseq", "Batch", matrix.su@meta.data$orig.ident)
matrix.su@meta.data$alignRate <- sampleTable[match(matrix.su@meta.data$orig.ident2, sampleTable$sampleLabels),]$uniq.align.rate
matrix.su@meta.data$Experiment <- sampleTable[match(matrix.su@meta.data$orig.ident2, sampleTable$sampleLabels),]$CellType
matrix.su@meta.data$Batch <- gsub( "_.*", "", matrix.su@meta.data$orig.ident)
View(matrix.su@meta.data)
unique(matrix.su@meta.data$Batch)
table(matrix.su@meta.data[,c(6,7)])
table(matrix.su@meta.data[,c(1,7)])

nrow(matrix.su@meta.data) #  308035 -4.9Gb
min(matrix.su@meta.data$nCount_RNA) # 259
min(matrix.su@meta.data$nFeature_RNA) # 240


matrix.su2 <- subset(matrix.su, subset = orig.ident2 %in% sampleTable_0.40$sampleLabels)
nrow(matrix.su2@meta.data) #  


keep_cells <- rownames(matrix.su@meta.data[matrix.su@meta.data$alignRate > 40 & matrix.su@meta.data$nFeature_RNA >= 300 ,])
length(keep_cells)  # 296897 with align rate > 40% (min 300g)
matrix.su <- SubsetData(matrix.su, cells = keep_cells)


nrow(matrix.su@meta.data) #  
table(matrix.su@meta.data[,c(6)])

keep_cells <- rownames(matrix.su@meta.data[matrix.su@meta.data$alignRate > 50 ,])
length(keep_cells)  #  294296 with align rate > 50% (min 300g)
matrix.su_50_300g <- subset(matrix.su, cells = keep_cells)

keep_cells <- rownames(matrix.su@meta.data[matrix.su@meta.data$nFeature_RNA >= 400 ,])
length(keep_cells)  #  cells with aligh rate > 40% & >= 400g
matrix.su_40_400g <- subset(matrix.su, cells = keep_cells)

keep_cells <- rownames(matrix.su_50_300g@meta.data[matrix.su_50_300g@meta.data$nFeature_RNA >= 400 ,])
length(keep_cells)  #  cells with aligh rate > 40% & >= 400g
matrix.su_50_400g <- subset(matrix.su, cells = keep_cells)

nrow(matrix.su_50_300g@meta.data) #  
table(matrix.su_50_300g@meta.data[,c(6)])

nrow(matrix.su_40_400g@meta.data) #  
table(matrix.su_40_400g@meta.data[,c(6)])

nrow(matrix.su_50_400g@meta.data) #  
table(matrix.su_50_400g@meta.data[,c(6)])

#saveRDS(matrix.su_50_300g, "CTR_dscj1_Realigned_50bp_matrix.su_alignRate50_300g.Rds") # cells
#saveRDS(matrix.su_40_400g, "CTR_dscj1_Realigned_50bp_matrix.su_alignRate40_400g.Rds") # cells
#saveRDS(matrix.su_50_400g, "CTR_dscj1_Realigned_50bp_matrix.su_alignRate50_400g.Rds") # cells







colnames(matrix.su) <- RenameCells(matrix.su, new.names = rownames(matrix.su@meta.data))


#write.csv(Sample_Cell_Genes.tbl, "Sample_Cell_Genes.tbl_172_samples.csv")


matrix.su@meta.data$orig.ident <- gsub("ForMatrix\\/", "", matrix.su@meta.data$orig.ident)
matrix.su@meta.data$orig.ident <- gsub("Dropseq_", "Dropseq", matrix.su@meta.data$orig.ident)
matrix.su@meta.data$orig.ident <- gsub("\\/SLX-", "_", matrix.su@meta.data$orig.ident)


rownames(matrix.su@meta.data) <- gsub("ForMatrix\\/", "", rownames(matrix.su@meta.data))
rownames(matrix.su@meta.data) <- gsub(".*\\/SLX-", "_", rownames(matrix.su@meta.data))
rownames(matrix.su@meta.data) <- gsub("_1_1_1_1_1_1", "", rownames(matrix.su@meta.data))


matrix.su@meta.data$cell_name <- rownames(matrix.su@meta.data)
matrix.su@meta.data$cell_name <- gsub(".*\\/SLX", "", matrix.su@meta.data$cell_name)
#matrix.su@meta.data$cell_name <- gsub("_1_1", "", matrix.su@meta.data$cell_name)
matrix.su@meta.data$cell_name <- gsub("-[0-9]{5}_N[0-9]{3}", "", matrix.su@meta.data$cell_name)
#matrix.su@meta.data$cell_name <- gsub("_[1-2]", "", matrix.su@meta.data$cell_name)
#matrix.su@meta.data$cell_name <- gsub("_", "", matrix.su@meta.data$cell_name)
matrix.su@meta.data$cell_name <- paste0(matrix.su@meta.data$orig.ident,"_", matrix.su@meta.data$cell_name)
rownames(matrix.su@meta.data) <- matrix.su@meta.data$cell_name
View(matrix.su@meta.data)
matrix.su@meta.data$cell_name <- gsub("_1_1_1_1", "", matrix.su@meta.data$cell_name)
matrix.su@meta.data$dups <- duplicated(matrix.su@meta.data$cell_name)


# change from seurat 2.3.4 to 3.0.2 to update the object:::
# 
#detach("package:Seurat", unload=TRUE)
#library("Seurat", lib.loc="/usr/local/lib/R/site-library")

#matrix.su <- UpdateSeuratObject(matrix.su)
#saveRDS(matrix.su, file = paste(resDir, "/matrix.su__181_samples_seurat.v.3.__400g.rds", sep = ""))

saveRDS(matrix.su, file = paste(resDir, "/matrix.su__172_samples_seurat.v.3.__300g.rds", sep = ""))




### Annotate new matrix.su :::


head(matrix.su@meta.data)
matrix.su@meta.data$Batch <- rownames(matrix.su@meta.data)
matrix.su@meta.data$Batch <- gsub("_.*", "", matrix.su@meta.data$Batch)
matrix.su@meta.data$Batch <- gsub("Batch", "Dropseq_", matrix.su@meta.data$Batch)
table(matrix.su@meta.data$Batch)



matrix.su@meta.data$Sample <- rownames(matrix.su@meta.data)
matrix.su@meta.data$Sample <- gsub("_[ACGTN]{12}", "", matrix.su@meta.data$Sample)
matrix.su@meta.data$BiolSample <- gsub("Batch[0-9]{4}_", "", matrix.su@meta.data$Sample)
length(unique(matrix.su@meta.data$BiolSample)) # 134 ???

matrix.su@meta.data$Experiment <- matrix.su@meta.data$Sample

samples_0.0 <- subset(sampleTable, sampleTable$CellType == "0.0")$sampleLabels
samples_I.1 <- subset(sampleTable, sampleTable$CellType == "I.1")$sampleLabels
samples_I.4 <- subset(sampleTable, sampleTable$CellType == "I.4")$sampleLabels
samples_I.24 <- subset(sampleTable, sampleTable$CellType == "I.24")$sampleLabels
samples_I.36 <- subset(sampleTable, sampleTable$CellType == "I.36")$sampleLabels
samples_I.48 <- subset(sampleTable, sampleTable$CellType == "I.48")$sampleLabels
samples_R.1 <- subset(sampleTable, sampleTable$CellType == "R.1")$sampleLabels
samples_R.4 <- subset(sampleTable, sampleTable$CellType == "R.4")$sampleLabels
samples_R.24 <- subset(sampleTable, sampleTable$CellType == "R.24")$sampleLabels
samples_R.36 <- subset(sampleTable, sampleTable$CellType == "R.36")$sampleLabels
samples_R.48 <- subset(sampleTable, sampleTable$CellType == "R.48")$sampleLabels

matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_0.0, "0.0", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_I.1, "I.1", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_I.4, "I.4", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_I.24, "I.24", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_I.36, "I.36", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_I.48, "I.48", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_R.1, "R.1", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_R.4, "R.4", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_R.24, "R.24", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_R.36, "R.36", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- ifelse(matrix.su@meta.data$orig.ident %in% samples_R.48, "R.48", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment <- gsub("15988_N704", "0.0", matrix.su@meta.data$Experiment)
matrix.su@meta.data$Experiment_batch <- paste("t",matrix.su@meta.data$Experiment, matrix.su@meta.data$Batch, sep = "_")

table(matrix.su@meta.data$Experiment)
names(table(matrix.su@meta.data$Batch))
head(matrix.su@meta.data)



metadata <- matrix.su@meta.data
keep_cells_batch0004rm <- rownames(metadata[metadata$Batch == "Dropseq_0014" | metadata$Batch == "Dropseq_0006" | metadata$Batch == "Dropseq_0007" | metadata$Batch == "Dropseq_0008" | metadata$Batch == "Dropseq_0011",])
length(keep_cells_batch0004rm)  # [1] 3360    247003
matrix.su <- subset(matrix.su, cells = keep_cells_batch0004rm)
length(unique(matrix.su@meta.data$orig.ident)) # 172
length(unique(matrix.su@meta.data$Sample)) # 172


#saveRDS(matrix.su, file = paste(resDir, "/matrix.su__172_orig.ident__seurat.v.3.__400g__0004rm.rds", sep = ""))
saveRDS(matrix.su, file = paste(resDir, "/matrix.su__172_orig.ident__seurat.v.3.__300g__0004rm.rds", sep = ""))








#matrix.su <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/SEURAT_172_samples_0006_0007_0008_0011_0014/matrix.su__172_orig.ident__seurat.v.3.__400g__0004rm.rds")
table(matrix.su@meta.data$Experiment[matrix.su@meta.data$Batch == "Dropseq_0006"])
table(matrix.su@meta.data$Experiment[matrix.su@meta.data$Batch == "Dropseq_0007"])
table(matrix.su@meta.data$Experiment[matrix.su@meta.data$Batch == "Dropseq_0008"])
table(matrix.su@meta.data$Experiment[matrix.su@meta.data$Batch == "Dropseq_0011"])
table(matrix.su@meta.data$Experiment[matrix.su@meta.data$Batch == "Dropseq_0014"])






matrix.su_backup <- matrix.su
matrix.su2 <- subset(x = matrix.su, subset = nFeature_RNA > 400)

cells_to_keep1 <- rownames(matrix.su2@meta.data[matrix.su2@meta.data$nFeature_RNA > 400,])
matrix.su_400 <- SubsetData(object = matrix.su, cells = cells_to_keep )
nrow(matrix.su_400@meta.data) # 60265
cells_to_keep400 <- cells_to_keep


cells_to_keep <- list()
for (i in 1:length(Experiment_batch_names)) {
  sample_expt_tmp <- subset(tmp_tbl, tmp_tbl$Experiment == Experiment_batch_names[i])
  sample_expt_tmp$cells <- rownames(sample_expt_tmp)
  if (nrow(sample_expt_tmp) > 10000){
    cells_to_keepp <- sample_n(sample_expt_tmp,  10000 ) #  as.integer(nrow(sample_expt_tmp)*0.8)
    cells_to_keep[[i]] <- cells_to_keepp$cells
  } else {  cells_to_keep[[i]] <- sample_expt_tmp$cells  }
}
cells_to_keep <- unlist(cells_to_keep)



matrix.su@meta.data$orig.ident2 <-matrix.su@meta.data$orig.ident
matrix.su@meta.data$orig.ident <-matrix.su@meta.data$Experiment
head(matrix.su@meta.data)
length(unique(matrix.su@meta.data$Sample))
length(unique(matrix.su@meta.data$Experiment))
(table(matrix.su@meta.data$Sample))
(table(matrix.su@meta.data$Experiment))



#cells_to_keep_0.5_list <- cells_to_keep
cells_to_keep_0.5_list <- cells_to_keep_0.5_list[c(1,2,4)]
#cells_to_keep_0.3_list <- cells_to_keep
cells_to_keep_0.3_list <- cells_to_keep_0.3_list[c(3,5)]
cells_to_keep_0.3_list <- unlist(cells_to_keep_0.3_list)
cells_to_keep_0.5_list <- unlist(cells_to_keep_0.5_list)

cells_to_keep <- unique(append(cells_to_keep_0.5_list, cells_to_keep_0.3_list))

cells_to_keep_0.6 <- unique(append(cells_to_keep400, cells_to_keep))
matrix.su_subsetted <- SubsetData(object = matrix.su, cells = cells_to_keep_0.6 )
nrow(matrix.su_subsetted@meta.data)


matrix.su_down <- subset(x = matrix.su, downsample = 7000)
cells_to_keep3 <- rownames(matrix.su_down@meta.data[matrix.su_down@meta.data$nFeature_RNA > 400,])
matrix.su_400x <- SubsetData(object = matrix.su_down, cells = cells_to_keep3 )


table(matrix.su_subsetted@meta.data$Experiment)
table(matrix.su@meta.data$Experiment)
table(matrix.su2@meta.data$Experiment)
table(matrix.su_400@meta.data$Experiment)
table(matrix.su_down@meta.data$Experiment)
table(matrix.su_down@meta.data$orig.ident2)
nrow(matrix.su_down@meta.data)
table(matrix.su_400x@meta.data$Experiment)
nrow(matrix.su_400x@meta.data)
table(matrix.su_400x@meta.data$orig.ident2)
min(matrix.su_down@meta.data$nCount_RNA)

saveRDS(matrix.su_down, file = paste(resDir, "/matrix.su_Dropseq_0011_R_15samples_min.gen50_94kcells.rds", sep = ""))

matrix.su_down <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/matrix.su_Dropseq_0011_R_15samples_min.gen50_94kcells.rds")  


