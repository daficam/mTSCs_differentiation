rm(list=ls())
gc()
options(bitmapType='cairo')


library("ggplot2")
library('VennDiagram');
library("data.table")
library("cowplot")
library("reshape")
library("ggrepel")
library("ggdendro")
library("useful")
library("MASS")
library("viridis")
library("scales")
library("biomaRt")
require("plyr")
library("ggalt")
library("ggridges")
library("dplyr")


baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)

hmisc_mat_r2 <- read.csv("r2_sample-to-sample_sc-to-bulk_165x165.csv")
hmisc_r2_molten <- reshape2::melt(hmisc_mat_r2[])
hmisc_r2_molten <- hmisc_r2_molten[hmisc_r2_molten$value != 1,]
hmisc_r2_molten <- hmisc_r2_molten[hmisc_r2_molten$value > 0.9,] # 4854
hmisc_r2_molten <- hmisc_r2_molten[order(hmisc_r2_molten$value, decreasing = T),]
hmisc_r2_molten <- hmisc_r2_molten[!duplicated(hmisc_r2_molten$value),] # 2427 sample pairs
colnames(hmisc_r2_molten) <- c("sample_1", "sample_2", "r2")

hmisc_mat_r2_x <- hmisc_mat_r2
rownames(hmisc_mat_r2_x) <- hmisc_mat_r2_x$X
hmisc_mat_r2_x <- hmisc_mat_r2_x[,-1]
colnames(hmisc_mat_r2_x) <- gsub( "_I.*", "", colnames(hmisc_mat_r2_x))
colnames(hmisc_mat_r2_x) <- gsub( "_R.*", "", colnames(hmisc_mat_r2_x))
colnames(hmisc_mat_r2_x) <- gsub( "_0.*", "", colnames(hmisc_mat_r2_x))
rownames(hmisc_mat_r2_x) <- gsub( "_I.*", "", rownames(hmisc_mat_r2_x))
rownames(hmisc_mat_r2_x) <- gsub( "_R.*", "", rownames(hmisc_mat_r2_x))
rownames(hmisc_mat_r2_x) <- gsub( "_0.*", "", rownames(hmisc_mat_r2_x))


uniq.align.rates <- read.csv("uniq.align.rates.csv", header = F)
colnames(uniq.align.rates) <- c("Batch","SLX", "SampleName", "AlignmentRate")
uniq.align.rates$Batch2 <- gsub( "CTR_DropSeq_", "Batch", uniq.align.rates$Batch)
uniq.align.rates <- subset(uniq.align.rates, uniq.align.rates$Batch2 %in% c( "Batch0004", "Batch0006", "Batch0007", "Batch0008" , "Batch0011", "Batch0014"))
uniq.align.rates$sampleLabels <- gsub( "SLX-", "", uniq.align.rates$SampleName)
uniq.align.rates$sampleLabels <- gsub( "\\.[ATXGC]{8}\\.", "_", uniq.align.rates$sampleLabels)
uniq.align.rates$sampleLabels <- paste(uniq.align.rates$Batch2 , uniq.align.rates$sampleLabels, sep= "_")

sampleTable <- read.csv("sampleTable_183samples_all_batches_alignrate.csv", row.names = 1)
#sampleTable <- read.csv("sampleTable_183samples_all_batches.csv")
#sampleTable$uniq.align.rate <- uniq.align.rates[match(sampleTable$sampleLabels, uniq.align.rates$sampleLabels),]$AlignmentRate
#sampleTable <- subset(sampleTable, sampleTable$uniq.align.rate > 40)
#sampleTable$fileUMI <- gsub("final.dge.txt.gz", "final.gatheredBarcodeDist.txt.gz", sampleTable$fileNameDGE)
FilesUMI <- sampleTable$fileUMI
names(FilesUMI) <- sampleTable$sampleLabels
#sampleTable[sampleTable$sampleLabels == "Batch0006_7635_N726",]$CellType <- "R.4"
#sampleTable[sampleTable$sampleLabels == "Batch0011_15980_N727",]$CellType <- "I.4"
#sampleTable[sampleTable$sampleLabels == "Batch0011_15987_N715",]$CellType <- "R.4"
#sampleTable$Treatment <- gsub( "\\..*", "", sampleTable$CellType)
#sampleTable$Age <- gsub( ".*\\.", "", sampleTable$CellType)
#write.csv(sampleTable, "sampleTable_183samples_all_batches_alignrate.csv")


# defile sample1 and sample2::

plot_shared_UMIs <- function(sample1, sample2) {
  
  file1 <- FilesUMI[names(FilesUMI) %in% sample1]
  file2 <- FilesUMI[names(FilesUMI) %in% sample2]
  
  mu.sample <- read.table(gzfile(file1), sep="\t",header=TRUE,  stringsAsFactors = TRUE) # row.names=1,
  mu.sample <- mu.sample[mu.sample$Gene == "Actb",]
  gt.sample <- read.table(gzfile(file2), sep="\t",header=TRUE,  stringsAsFactors = TRUE) # row.names=1,
  gt.sample <- gt.sample[gt.sample$Gene == "Actb",]
  colnames(mu.sample) <- c("UMI", "Gene", "CELLBARCODE", "COUNT")
  colnames(gt.sample) <- c("UMI", "Gene", "CELLBARCODE", "COUNT")
  
  gt.sample.agg           <- aggregate(gt.sample$COUNT, by=list(UMI=gt.sample$UMI), FUN=sum)
  mu.sample.agg           <- aggregate(mu.sample$COUNT, by=list(UMI=mu.sample$UMI), FUN=sum)
  colnames(gt.sample.agg) <- c("UMI", "gt.sample")
  colnames(mu.sample.agg) <- c("UMI", "mu.sample")
  corr.df                 <- merge(gt.sample.agg, mu.sample.agg, by="UMI")
  n_common_UMIs           <- nrow(corr.df)
  
  plt.corr <- ggplot(corr.df, aes(x = gt.sample, y=mu.sample, col= 1-(abs(mu.sample - gt.sample)/(mu.sample + gt.sample)) ) ) +
    geom_abline(slope=1, intercept=0) +
    geom_point(size = 2, alpha=0.5 ) + 
    scale_color_continuous(name="Correlation") +
    ggtitle("UMI overlap") +
    theme(aspect.ratio=1) + 
    labs(y=sample2, x = sample1)
  return(plt.corr) # 
}


calc_n_shared_UMIs <- function(sample1, sample2) {
  
  file1 <- FilesUMI[names(FilesUMI) %in% sample1]
  file2 <- FilesUMI[names(FilesUMI) %in% sample2]
  
  mu.sample <- read.table(gzfile(file1), sep="\t",header=TRUE,  stringsAsFactors = TRUE) # row.names=1,
  mu.sample <- mu.sample[mu.sample$Gene == "Actb",]
  gt.sample <- read.table(gzfile(file2), sep="\t",header=TRUE,  stringsAsFactors = TRUE) # row.names=1,
  gt.sample <- gt.sample[gt.sample$Gene == "Actb",]
  colnames(mu.sample) <- c("UMI", "Gene", "CELLBARCODE", "COUNT")
  colnames(gt.sample) <- c("UMI", "Gene", "CELLBARCODE", "COUNT")
  
  gt.sample.agg           <- aggregate(gt.sample$COUNT, by=list(UMI=gt.sample$UMI), FUN=sum)
  mu.sample.agg           <- aggregate(mu.sample$COUNT, by=list(UMI=mu.sample$UMI), FUN=sum)
  colnames(gt.sample.agg) <- c("UMI", "gt.sample")
  colnames(mu.sample.agg) <- c("UMI", "mu.sample")
  corr.df                 <- merge(gt.sample.agg, mu.sample.agg, by="UMI")
  n_common_UMIs           <- nrow(corr.df)
  
  return(n_common_UMIs) # plt.corr
}


# Mixed up files (mu):::
#FilesUMI_mu <- FilesUMI[names(FilesUMI) %in% c("Batch0011_15980_N727", "Batch0011_15987_N715", "Batch0011_15988_N704")]

tmp_df <- as.data.frame((hmisc_mat_r2[hmisc_mat_r2$X %in% c("Batch0011_15980_N727_15980_xN727", "Batch0011_15987_N715_15987_xN715", "Batch0011_15988_N704_0.0"),]))
rownames(tmp_df) <- tmp_df$X
tmp_df <- tmp_df[,-1]
tmp_df <- as.data.frame(t(tmp_df))

# example of true replicates:
plt.corr_N704  <- plot_shared_UMIs("Batch0011_15988_N704",  "Batch0006_7633_N704")
plt.corr_N715  <- plot_shared_UMIs("Batch0011_15987_N715",  "Batch0006_7633_N715" )
plt.corr_N727  <- plot_shared_UMIs("Batch0011_15980_N727",  "Batch0006_7635_N727")

plot_grid(plt.corr_N704, plt.corr_N715, plt.corr_N727, ncol=3, vhjust =1)

# example of non-replicates:
plt.corr_N704x  <- plot_shared_UMIs("Batch0011_15988_N704",  "Batch0014_9616_N723")
plt.corr_N715x  <- plot_shared_UMIs("Batch0011_15987_N715",  "Batch0011_15988_N716" )
plt.corr_N727x  <- plot_shared_UMIs("Batch0011_15980_N727",  "Batch0006_7635_N705")

plot_grid(plt.corr_N704x, plt.corr_N715x, plt.corr_N727x, ncol=3, vhjust =1)



find_all_pairs <- function(sample_list){
  x_pairs <- purrr::cross2(seq_along(1:length(sample_list)), seq_along(1:length(sample_list)))
  x_pairs <- as.data.frame(t(matrix(unlist(x_pairs), nrow = 2)))
  x_pairs <- x_pairs[!x_pairs$V1 == x_pairs$V2 ,]
  x_pairs <- x_pairs[order(x_pairs$V1),]
  rownames(x_pairs) <- paste0("samples_", x_pairs$V1, "_", x_pairs$V2)
  x_pairs <- as.data.frame(t(x_pairs))
  return(x_pairs)
}


#sample_list <- Dups_setA
  
find_replicates_among_samples <- function(sample_list){
  x_pairs <- find_all_pairs(sample_list)
  plt.corr_x  <- calc_n_shared_UMIs(sample_list[x_pairs[,1][[1]]],  sample_list[x_pairs[,1][[2]]])
  
  UMI_n_x_list <- list()
  for(i in seq_along(colnames(x_pairs))){
    plt.corr_x  <- calc_n_shared_UMIs(sample_list[x_pairs[,i][[1]]],  sample_list[x_pairs[,i][[2]]])
    UMI_n_x_list[[i]]    <- plt.corr_x
  }
  
  shared_UMI_tbl_pairs <- as.data.frame(t(x_pairs))
  shared_UMI_tbl_pairs$no_shared_UMIs <- unlist(UMI_n_x_list)
  names(sample_list) <- seq_along(sample_list)
  shared_UMI_tbl_pairs$sample1 <- sample_list[match(shared_UMI_tbl_pairs$V1, names(sample_list))]
  shared_UMI_tbl_pairs$sample2 <- sample_list[match(shared_UMI_tbl_pairs$V2, names(sample_list))]
  shared_UMI_tbl_pairs <- shared_UMI_tbl_pairs[order(shared_UMI_tbl_pairs$no_shared_UMIs, decreasing = T),]
  shared_UMI_tbl_pairs <- shared_UMI_tbl_pairs[shared_UMI_tbl_pairs$no_shared_UMIs > 50,]
  return(shared_UMI_tbl_pairs[,-c(1,2)])
}



Dups_setA <- c("Batch0004_9617_N720", "Batch0004_9617_N721", "Batch0004_9617_N724", "Batch0004_9617_N726", "Batch0006_7633_N719","Batch0006_7635_N719", "Batch0007_7637_N724")
Dups_setB <- c("Batch0004_9617_N722", "Batch0004_9617_N719", "Batch0004_9617_N718", "Batch0004_9617_N727","Batch0004_9617_N729", "Batch0006_7633_N704","Batch0006_7635_N704")
Dups_1A <- c("Batch0006_7633_N723", "Batch0006_7635_N723", "Batch0008_9339_N723", "Batch0008_9616_N723", "Batch0014_9339_N723","Batch0014_9616_N723")
Dups_2B <- c("Batch0006_7633_N712", "Batch0006_7635_N712", "Batch0008_9339_N712", "Batch0008_9616_N712", "Batch0011_15988_N712","Batch0014_9339_N712", "Batch0014_9616_N712")
Dups_2C <- c("Batch0006_7633_N703", "Batch0006_7635_N703", "Batch0011_15985_N703")
Dups_2D <- c("Batch0006_7633_N718", "Batch0006_7635_N718", "Batch0007_7637_N718", "Batch0008_9339_N718", "Batch0008_9616_N718","Batch0011_15988_N718", "Batch0014_9339_N718", "Batch0014_9616_N718")
Dups_1B <- c("Batch0006_7633_N711", "Batch0006_7635_N711", "Batch0007_7637_N711", "Batch0008_9339_N711", "Batch0008_9616_N711","Batch0011_15987_N711", "Batch0014_9339_N711", "Batch0014_9616_N711")
Dups_1C <- c("Batch0006_7633_N724", "Batch0006_7635_N724", "Batch0008_9339_N724", "Batch0008_9616_N724","Batch0011_15980_N724", "Batch0014_9339_N724", "Batch0014_9616_N724")
Dups_1D <- c("Batch0006_7633_N716", "Batch0006_7635_N716", "Batch0008_9339_N716", "Batch0008_9616_N716","Batch0011_15988_N716", "Batch0014_9339_N716", "Batch0014_9616_N716")
Dups_5B <- c("Batch0008_9616_N714", "Batch0014_9339_N714", "Batch0014_9616_N714")
Dups_5C <- c("Batch0006_7633_N705", "Batch0006_7635_N705", "Batch0011_15987_N705")
Dups_6A <- c("Batch0006_7635_N727", "Batch0008_9339_N727", "Batch0008_9616_N727", "Batch0011_15980_N727", "Batch0014_9339_N727","Batch0014_9616_N727")
Dups_5A <- c("Batch0006_7633_N726", "Batch0008_9339_N726", "Batch0008_9616_N726", "Batch0011_15980_N726", "Batch0014_9339_N726","Batch0014_9616_N726", "Batch0006_7635_N726")
Dups_6B <- c("Batch0006_7633_N715", "Batch0006_7635_N715", "Batch0008_9339_N715", "Batch0008_9616_N715", "Batch0014_9339_N715","Batch0014_9616_N715", "Batch0011_15987_N715")
Dups_6C <- c("Batch0006_7633_N720", "Batch0006_7635_N720", "Batch0011_15982_N720")
Dups_6D <- c("Batch0008_9616_N720", "Batch0014_9616_N720")
Dups_7A <- c("Batch0006_7633_N728", "Batch0006_7635_N728", "Batch0007_7637_N728", "Batch0008_9339_N728", "Batch0008_9616_N728","Batch0011_15982_N728", "Batch0014_9339_N728", "Batch0014_9616_N728")
Dups_8B <- c("Batch0006_7633_N710", "Batch0006_7635_N710", "Batch0008_9339_N710", "Batch0008_9616_N710", "Batch0011_15988_N710","Batch0014_9339_N710", "Batch0014_9616_N710")
Dups_8C <- c("Batch0006_7633_N706", "Batch0006_7635_N706", "Batch0011_15987_N706")
Dups_8D <- c("Batch0006_7633_N722", "Batch0006_7635_N722", "Batch0008_9339_N722", "Batch0008_9616_N722", "Batch0011_15980_N722","Batch0014_9339_N722", "Batch0014_9616_N722")
Dups_7B <- c("Batch0006_7633_N707", "Batch0006_7635_N707", "Batch0008_9339_N707", "Batch0008_9616_N707", "Batch0011_15980_N707","Batch0014_9339_N707", "Batch0014_9616_N707")
Dups_8A <- c("Batch0006_7633_N729", "Batch0006_7635_N729", "Batch0008_9339_N729", "Batch0008_9616_N729", "Batch0011_15980_N729","Batch0014_9339_N729", "Batch0014_9616_N729")
Dups_36I1 <- c("Batch0007_7637_N701", "Batch0011_15985_N701")
Dups_36I2 <- c("Batch0007_7637_N702", "Batch0011_15987_N702")
Dups_36I3 <- c("Batch0007_7637_N703", "Batch0011_15987_N703")
Dups_36I4 <- c("Batch0007_7637_N704", "Batch0011_15985_N704")
Dups_36R1 <- c("Batch0007_7637_N706", "Batch0011_15985_N706")
Dups_36R2 <- c("Batch0007_7637_N707", "Batch0011_15985_N707")
Dups_36R3 <- c("Batch0007_7637_N710", "Batch0011_15987_N710")
Dups_36R4 <- c("Batch0007_7637_N712", "Batch0011_15982_N712")
Dups_48I1 <- c("Batch0007_7637_N705", "Batch0011_15985_N705")
Dups_48I2 <- c("Batch0007_7637_N716", "Batch0011_15982_N716")
Dups_48I3 <- c("Batch0007_7637_N719", "Batch0011_15988_N719")
Dups_48I4 <- c("Batch0007_7637_N720", "Batch0011_15988_N720")
Dups_48R1 <- c("Batch0007_7637_N729", "Batch0011_15985_N729")
Dups_48R2 <- c("Batch0007_7637_N722", "Batch0011_15982_N722")
Dups_48R3 <- c("Batch0007_7637_N723", "Batch0011_15980_N723")
Dups_48R4 <- c("Batch0007_7637_N726", "Batch0011_15982_N726")



Dups_setA_reps <- find_replicates_among_samples(Dups_setA)
Dups_setB_reps <- find_replicates_among_samples(Dups_setB)
Dups_1A_reps   <- find_replicates_among_samples(Dups_1A)
Dups_2B_reps   <- find_replicates_among_samples(Dups_2B)
Dups_2C_reps   <- find_replicates_among_samples(Dups_2C)
Dups_2D_reps   <- find_replicates_among_samples(Dups_2D)
Dups_1B_reps   <- find_replicates_among_samples(Dups_1B)
Dups_1C_reps   <- find_replicates_among_samples(Dups_1C)
Dups_1D_reps   <- find_replicates_among_samples(Dups_1D)
Dups_5B_reps   <- find_replicates_among_samples(Dups_5B)
Dups_5C_reps   <- find_replicates_among_samples(Dups_5C)
Dups_6A_reps   <- find_replicates_among_samples(Dups_6A)
Dups_5A_reps   <- find_replicates_among_samples(Dups_5A)
Dups_6B_reps   <- find_replicates_among_samples(Dups_6B)
Dups_6C_reps   <- find_replicates_among_samples(Dups_6C)
Dups_6D_reps   <- find_replicates_among_samples(Dups_6D)
Dups_7A_reps   <- find_replicates_among_samples(Dups_7A)
Dups_8B_reps   <- find_replicates_among_samples(Dups_8B)
Dups_8C_reps   <- find_replicates_among_samples(Dups_8C)
Dups_8D_reps   <- find_replicates_among_samples(Dups_8D)
Dups_7B_reps   <- find_replicates_among_samples(Dups_7B)
Dups_8A_reps   <- find_replicates_among_samples(Dups_8A)
Dups_36I1_reps <- find_replicates_among_samples(Dups_36I1)
Dups_36I2_reps <- find_replicates_among_samples(Dups_36I2)
Dups_36I3_reps <- find_replicates_among_samples(Dups_36I3)
Dups_36I4_reps <- find_replicates_among_samples(Dups_36I4)
Dups_36R1_reps <- find_replicates_among_samples(Dups_36R1)
Dups_36R2_reps <- find_replicates_among_samples(Dups_36R2)
Dups_36R3_reps <- find_replicates_among_samples(Dups_36R3)
Dups_36R4_reps <- find_replicates_among_samples(Dups_36R4)
Dups_48I1_reps <- find_replicates_among_samples(Dups_48I1)
Dups_48I2_reps <- find_replicates_among_samples(Dups_48I2)
Dups_48I3_reps <- find_replicates_among_samples(Dups_48I3)
Dups_48I4_reps <- find_replicates_among_samples(Dups_48I4)
Dups_48R1_reps <- find_replicates_among_samples(Dups_48R1)
Dups_48R2_reps <- find_replicates_among_samples(Dups_48R2)
Dups_48R3_reps <- find_replicates_among_samples(Dups_48R3)
Dups_48R4_reps <- find_replicates_among_samples(Dups_48R4)


plot_shared_UMIs(Dups_1A[1], Dups_1A[3])
plot_shared_UMIs(Dups_1A[1], Dups_1A[2])
plot_shared_UMIs(Dups_1A[4], Dups_1A[3])
plot_shared_UMIs(Dups_1A[2], Dups_1A[5])


#min(hmisc_mat_r2_x[colnames(hmisc_mat_r2_x) %in% Dups_2C_reps$sample1, rownames(hmisc_mat_r2_x) %in% Dups_2C_reps$sample2])

list_dup_reps <- list(Dups_setA_reps=Dups_setA_reps, Dups_setB_reps=Dups_setB_reps, Dups_1A_reps=Dups_1A_reps, 
                      Dups_2B_reps=Dups_2B_reps,Dups_2C_reps=Dups_2C_reps, Dups_2D_reps=Dups_2D_reps, 
                      Dups_1B_reps=Dups_1B_reps, Dups_1C_reps=Dups_1C_reps, Dups_1D_reps=Dups_1D_reps, 
                      Dups_5B_reps=Dups_5B_reps, Dups_5C_reps=Dups_5C_reps, Dups_6A_reps=Dups_6A_reps, 
                      Dups_5A_reps=Dups_5A_reps, Dups_6B_reps=Dups_6B_reps, Dups_6C_reps=Dups_6C_reps,
                      Dups_6D_reps=Dups_6D_reps, Dups_7A_reps=Dups_7A_reps, Dups_8B_reps=Dups_8B_reps, 
                      Dups_8C_reps=Dups_8C_reps, Dups_8D_reps=Dups_8D_reps, Dups_7B_reps=Dups_7B_reps, 
                      Dups_8A_reps=Dups_8A_reps, Dups_36I1_reps=Dups_36I1_reps, Dups_36I2_reps=Dups_36I2_reps,
                      Dups_36I3_reps=Dups_36I3_reps, Dups_36I4_reps=Dups_36I4_reps, Dups_36R1_reps=Dups_36R1_reps,
                      Dups_36R2_reps=Dups_36R2_reps, Dups_36R3_reps=Dups_36R3_reps, Dups_36R4_reps=Dups_36R4_reps, 
                      Dups_48I1_reps=Dups_48I1_reps, Dups_48I2_reps=Dups_48I2_reps, Dups_48I3_reps=Dups_48I3_reps,
                      Dups_48I4_reps=Dups_48I4_reps, Dups_48R1_reps=Dups_48R1_reps, Dups_48R2_reps=Dups_48R2_reps,
                      Dups_48R3_reps=Dups_48R3_reps, Dups_48R4_reps=Dups_48R4_reps )

mn_r2 <- list()
for(i in seq_along(list_dup_reps)) {
  df <- list_dup_reps[[i]]
  mn_r2[[i]] <- min(hmisc_mat_r2_x[colnames(hmisc_mat_r2_x) %in% df$sample1, rownames(hmisc_mat_r2_x) %in% df$sample2])
}
min(unlist(mn_r2))

sure_reps <- data.frame(r2 = unlist(mn_r2))
rownames(sure_reps) <- names(list_dup_reps)


### ASSESS unsure reps:::

unsure_reps <- subset(sure_reps, sure_reps$r2 < 0.85)
sure_reps   <- subset(sure_reps, sure_reps$r2 >= 0.85)
rownames(unsure_reps) # "Dups_2B_reps" "Dups_1B_reps" "Dups_1C_reps"

Dups_2B_reps <- read.csv("Dups_reps/Dups_2B_reps.csv", row.names = 1)
Dups_1B_reps <- read.csv("Dups_reps/Dups_1B_reps.csv", row.names = 1)
Dups_1C_reps <- read.csv("Dups_reps/Dups_1C_reps.csv", row.names = 1)

mat1 <- hmisc_mat_r2_x[colnames(hmisc_mat_r2_x) %in% Dups_2B_reps$sample1, rownames(hmisc_mat_r2_x) %in% Dups_2B_reps$sample2]
mat2 <- hmisc_mat_r2_x[colnames(hmisc_mat_r2_x) %in% Dups_1B_reps$sample1, rownames(hmisc_mat_r2_x) %in% Dups_1B_reps$sample2]
mat3 <- hmisc_mat_r2_x[colnames(hmisc_mat_r2_x) %in% Dups_1C_reps$sample1, rownames(hmisc_mat_r2_x) %in% Dups_1C_reps$sample2]
mat4_sure <- hmisc_mat_r2_x[colnames(hmisc_mat_r2_x) %in% Dups_7B_reps$sample1, rownames(hmisc_mat_r2_x) %in% Dups_7B_reps$sample2]


#df_unsure1_r2_molten <- reshape2::melt(hmisc_mat_r2_x[colnames(hmisc_mat_r2_x) %in% Dups_2B_reps$sample1, rownames(hmisc_mat_r2_x) %in% Dups_2B_reps$sample2])

f2 = colorRamp2( c(0, 0.8499,0.85, 0.95, 1), c("white", "white","lightskyblue",  "deepskyblue3", "blue4"), space = "RGB") 
ht3 = Heatmap(as.matrix(mat1),  col = f2, name = "sc-to-Bulk",  row_title = "", column_title = "sc-to-Bulk duplicates", show_row_names = TRUE, heatmap_legend_param = list(title = "Correlation", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = T, cluster_rows = T , row_names_side ="left")
ht3



tbl_for_merging <- data.frame(no_shared_UMIs=numeric(0), sample1=character(0), sample2=character(0), rep_cluster=character(0))

list_dup_reps_pairs <- list()
for ( i in seq_along(list_dup_reps)) {
  df <- list_dup_reps[[i]]
  df <- df[order(df$no_shared_UMIs),]
  df <- df[!duplicated(df$no_shared_UMIs),]
  df$rep_clust <- names(list_dup_reps)[[i]]
  list_dup_reps_pairs[[i]] <- df
  tbl_for_merging <- rbind(tbl_for_merging, df)
}

tbl_for_merging$sampleFile1 <- sampleTable[match(tbl_for_merging$sample1, sampleTable$sampleLabels),]$fileNameDGE
tbl_for_merging$sampleFile2 <- sampleTable[match(tbl_for_merging$sample2, sampleTable$sampleLabels),]$fileNameDGE
head(tbl_for_merging,2)




tbl_for_merging[tbl_for_merging$rep_clust %in% rownames(unsure_reps),]


hmisc_mat_molten <- reshape2::melt(hmisc_mat_r2)
colnames(hmisc_mat_molten) <- c("sample1", "sample2", "r2")
hmisc_mat_molten$sample1 <- gsub("_0.*", "", hmisc_mat_molten$sample1)
hmisc_mat_molten$sample1 <- gsub("_R.*", "", hmisc_mat_molten$sample1)
hmisc_mat_molten$sample1 <- gsub("_I.*", "", hmisc_mat_molten$sample1)
hmisc_mat_molten$sample2 <- gsub("_0.*", "", hmisc_mat_molten$sample2)
hmisc_mat_molten$sample2 <- gsub("_R.*", "", hmisc_mat_molten$sample2)
hmisc_mat_molten$sample2 <- gsub("_I.*", "", hmisc_mat_molten$sample2)

hmisc_mat_molten$match_name <- paste(hmisc_mat_molten$sample1, hmisc_mat_molten$sample2, sep = "_")
tbl_for_merging$match_name <- paste(tbl_for_merging$sample1, tbl_for_merging$sample2, sep = "_")
tbl_for_merging$r2 <- hmisc_mat_molten[match(tbl_for_merging$match_name, hmisc_mat_molten$match_name),]$r2

tbl_for_merging <- subset(tbl_for_merging, tbl_for_merging$r2 > 0.85)
#write.csv(tbl_for_merging, "tbl_for_merging.csv")
tbl_for_merging <- read.csv("tbl_for_merging_with_unsure.csv")

matrix.su <- readRDS("CTR_dscj1_Realigned_50bp_matrix.su_alignRate40_300g_added2unknownSamples.Rds")

head(matrix.su@meta.data,2)
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0006_7635_N724",]) # 1622 cells
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0006_7633_N711",]) # 1660 cells
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0008_9339_N711",]) # 333 cells
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident2 == "Batch0014_9339_N711",]) # 366 cells

Dups_check_reps <- find_replicates_among_samples(c("Batch0014_9616_N712", "Batch0006_7635_N712"))




##### create lists of replicates to merge::::

head(tbl_for_merging,3)

list_of_reps_to_merge <- split(tbl_for_merging, tbl_for_merging$rep_clust)
names(list_of_reps_to_merge) <- unique(tbl_for_merging$rep_clust)
 
#list_reps <- list()
#list_reps_files <- list() 


#reps <- list_of_reps_to_merge[[1]][1,c(3,4)]

#for (i in seq_along(list_of_reps_to_merge)){
#  list_reps[[i]] <- unique(c(as.character(list_of_reps_to_merge[[i]]$sample1), as.character(list_of_reps_to_merge[[i]]$sample2)))
#  list_reps_files[[i]] <- unique(c(as.character(list_of_reps_to_merge[[i]]$sampleFile1), as.character(list_of_reps_to_merge[[i]]$sampleFile2)))
#}

names(list_reps) <- unique(tbl_for_merging$rep_clust)
names(list_reps_files) <- unique(tbl_for_merging$rep_clust)

saveRDS(list_reps, "list_reps.Rds")
saveRDS(list_reps_files, "list_reps_files.Rds")


final_set_to_merge <- read.csv("Final_set_to_merge.csv")
rownames(final_set_to_merge) <- final_set_to_merge$Name_of_rep_set
final_set_to_merge <- final_set_to_merge[,-c(10)]
#final_set_to_merge_t <- as.data.frame(t(final_set_to_merge))

final_set_to_merge$Name_of_rep_set <- as.character(final_set_to_merge$Name_of_rep_set)
final_set_to_merge <- final_set_to_merge[order(final_set_to_merge$Name_of_rep_set),]
names_sets <- unique(final_set_to_merge$Name_of_rep_set)
#final_set_to_merge$set <- final_set_to_merge$Name_of_rep_set
final_set_to_merge_list <- split(final_set_to_merge, final_set_to_merge$Name_of_rep_set)

#for(i in seq_along(final_set_to_merge_list)){
#  names(final_set_to_merge_list)[[i]] <- final_set_to_merge_list[[i]][,10]
#}


list_reps <- list()
list_reps_files <- list() 

for (i in seq_along(final_set_to_merge_list)){
  reps <- as.character(t(final_set_to_merge_list[[i]][,-c(1,11)]))
  #reps <- reps[reps != ""]
  list_reps[[i]] <- reps
  list_reps_files[[i]] <- sampleTable[match(reps, sampleTable$sampleLabels),]$fileNameDGE
}
names(list_reps) <- names_sets
names(list_reps_files) <- names_sets

list_reps2 <- list_reps
#saveRDS(list_reps, "Final_list_reps.Rds")
#saveRDS(list_reps_files, "Final_list_reps_files.Rds")
df1 <- data.frame(sample_name = list_reps)
df2 <- data.frame(files= list_reps_files)
