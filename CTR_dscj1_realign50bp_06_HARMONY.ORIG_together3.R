#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()
options(bitmapType='cairo')
#httr::set_config(httr::config(ssl_verifypeer = FALSE))


library(cowplot)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(clustree)
#library(harmony)
library(viridis)
library(colorspace)
library(optrees)
library(deldir)
#library(alphahull)

#options(future.globals.maxSize = 4000 * 1024^2)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
Project <- "CTR_dscj1_HARMONY_orig_500g"

reso <- "SCT_snn_res.0.8"

dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")

head(dropseq.integrated@meta.data,2)
Idents(dropseq.integrated) <- dropseq.integrated@meta.data[[reso]]


#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_clust24_25rm.Rds")
nrow(dropseq.integrated@meta.data[dropseq.integrated@meta.data$orig.ident2 == "Batch0004_9617_N721",])
nrow(dropseq.integrated@meta.data[dropseq.integrated@meta.data$orig.ident2 == "Batch0004_9617_N729",])


#saveRDS(matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/matrix.su_400g_SCT.Rds")

nrow(dropseq.integrated@meta.data)
mean(dropseq.integrated@meta.data$nCount_RNA)
mean(dropseq.integrated@meta.data$nFeature_RNA)




# FeaturePlot: Phf8, Gcm1, E2f8, ...
pdf(paste( "FeaturePlot", dataset, "_green_Phf8", ".pdf", sep="_"), width=12,height=8)
par(bg=NA)
FeaturePlot(dropseq.integrated,  features =c("Phf8","Gcm1","E2f8") )
dev.off()


tmp_tbl <- as.data.frame(table(dropseq.integrated@meta.data[,c(reso, "Experiment")]))
write.csv(tmp_tbl, "Harmony.orig_500g_together_TABLE_no_of_cells_per_Cluster_Experiment.csv")

colnames(dropseq.integrated@meta.data)
cell_table <- dropseq.integrated@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "Experiment",  "Age" , "Treatment",  "Batch"   )]
write.csv(cell_table, "Harmony.orig_500g_together_cell_table.csv")




message("+-------------------------------------------------------------------------------")
message("+               load in dataset for harmony integration                         ")
message("+-------------------------------------------------------------------------------")


#matrix.su <- readRDS("CTR_dscj1_Realigned_MERGED_50bp_matrix.su_300g_40pc_alignRate_53_SCTransformClust_UMAP_30L.Rds")
#matrix.su <- readRDS("HARMONY.orig_400g_500g/matrix.su_400g_SCT.Rds") # tsting if 400g or 500g cutoff works better???
#matrix.su <- readRDS("HARMONY.orig_400g_500g/matrix.su_500g_SCT.Rds") # tsting if 400g or 500g cutoff works better???


nrow(matrix.su@meta.data)             # 229793 
min(matrix.su@meta.data$nFeature_RNA) # 300
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident == "Batch0004_9617_N721",])
nrow(matrix.su@meta.data[matrix.su@meta.data$orig.ident == "Batch0004_9617_N729",])
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N721")
matrix.su <- subset(matrix.su, orig.ident != "Batch0004_9617_N729")

head(matrix.su@meta.data,2)
tmp_df <- unique(matrix.su@meta.data[,c("Experiment", "orig.ident2", "Batch")])
table(tmp_df[,c(1,3)])

matrix.su <- subset(matrix.su,  nFeature_RNA > 400) # 169094 cells
min(matrix.su@meta.data$nFeature_RNA)


message("+-------------------------------------------------------------------------------")
message("+                                HARMONY orig                               ")
message("+-------------------------------------------------------------------------------")

harmony_matrix.su <- RunHarmony(matrix.su, group.by.vars = "Batch", reduction = "pca", assay.use="SCT", reduction.save = "harmony") 
gc()
harmony_matrix.su <- FindNeighbors(harmony_matrix.su, reduction = "harmony", dims = 1:30) 
gc()
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.2)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.6)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 0.8)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 1.5)
harmony_matrix.su <- FindClusters(harmony_matrix.su, resolution = 2)

harmony_matrix.su <- RunUMAP(harmony_matrix.su, reduction = "harmony", dims = 1:30, min.dist = 0.01, n.neighbors = 25, repulsion.strength = 2, spread = 2L )

#saveRDS(harmony_matrix.su, "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")



message("+-------------------------------------------------------------------------------")
message("+                    load in  HARMONY object                                    ")
message("+-------------------------------------------------------------------------------")

#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729_clust24_25rm.Rds")
#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I/harmony_matrix.su_umap_n.neigh25_repuls2_spread2L_rm0004_9617_N729.Rds")
#dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_400g.Rds")
dropseq.integrated <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")

Project <- "CTR_dscj1_HARMONY_orig_500g"
gc()

# quick checks:

head(dropseq.integrated@meta.data,2)

dropseq.integrated@meta.data$Experiment <- factor(dropseq.integrated@meta.data$Experiment, levels = c("0.0","I.1","R.1","I.4","R.4","I.24","R.24","I.36","R.36","I.48","R.48"))
table(dropseq.integrated@meta.data$Experiment)

min(dropseq.integrated@meta.data$nFeature_RNA) # 500g

nrow(dropseq.integrated@meta.data) # 136291



message("+-------------------------------------------------------------------------------")
message("+              Explore and QC integrated dataset                                ")
message("+-------------------------------------------------------------------------------")

plots <- DimPlot(dropseq.integrated, group.by = c("Batch", "Experiment"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))


plots2 <- DimPlot(dropseq.integrated, group.by = "Batch", split.by = "Experiment", ncol = 3)


message("+-------------------------------------------------------------------------------")
message("+                   violin plots                   ")
message("+-------------------------------------------------------------------------------")

head(dropseq.integrated@meta.data,2)

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Batch
VlnPlot(dropseq.integrated, features = c("nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 1, pt.size = 0)

pdf(paste(Project, "harmony_500g___QC_VlnPlot_nCountRNA_nFeatRNA_percent.mt.pdf", sep=""), width=6,height=9)
par(bg=NA)
VlnPlot(dropseq.integrated, features = c("nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 1, pt.size = 0)
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                       CLUSTREE                                                ")
message("+-------------------------------------------------------------------------------")

head(dropseq.integrated@meta.data,2)
#dropseq.integrated@meta.data <- dropseq.integrated@meta.data[,-c(18,19)]

tree_to_plot <- clustree::clustree(dropseq.integrated, prefix = "SCT_snn_res.") # "RNA_snn_res."

pdf(paste(Project, "harmony_500g___matrix.su_clustree_0_to_2.pdf", sep=""), width=12,height=8)
par(bg=NA)
tree_to_plot
dev.off()



message("+-------------------------------------------------------------------------------")
message("+                       pengminshi/MRtree                                       ")
message("+-------------------------------------------------------------------------------")
# https://htmlpreview.github.io/?https://github.com/pengminshi/MRtree/blob/master/vignettes/MRtree-tutorial.html
# https://github.com/pengminshi/MRtree

library(mrtree)

dat <- dropseq.integrated

set.seed(42)

metadata = dat$metadata
rownames(metadata) = dat$metadata$cellid
ref.labels = dat$metadata$type

# specify the resolution parameters
# resolutions = seq(0.1, sqrt(3), 0.05)^2

# alternatively and preferrably, we provide a sampling tool to sample resolution parameters to uniformly cover different scales
A = seurat_get_nn_graph(counts=dat$counts, metadata=metadata, npc=10)
resolutions = modularity_event_sampling(A=A, n.res=30, gamma.min=0.01, gamma.max=2.5) # sample based on the similarity matrix

# clustering using Suerat 
seurat.out = sc_clustering.seurat(counts=dat$counts, resolutions=resolutions, metadata=metadata, npcs=10,
                                  min.cells=0, min.features=0, scale.factor=10000, return.seurat.object=T,
                                  vars.to.regress=NULL, find.variable.features=F, verbose=F)

# initial cluster tree from Seurat flat clustering
plot_clustree(labelmat=seurat.out$seurat.clusters, prefix ='SCT_snn_res.', 
              ref.labels = ref.labels, plot.ref = F)


#Then we apply MRtree to ubtain the hierarchical cluster tree, visualized using a dendrogram, with a pie chart on each tree node detailing the cluster composition given the known true labels.

out = mrtree(seurat.out$obj, consensus=F, augment.path=F)
# if there are few partitions per k, within resolution consensus step can speed up the algorithm
# weight per sample is encoraged if the classes are imbalanced

plot_tree(labelmat=out$labelmat.mrtree, ref.labels=ref.labels, plot.piechart = T,
          node.size = 0.4, tip.label.dist = 10, bottom.margin=30 )



# We evaluate te per-resolution clustering performance with a novel index adapted from Adjusted Rand Index to accrount for te bias for resolution.

ks.flat = apply(out$labelmat.flat, 2, FUN=function(x) length(unique(x)))
ks.mrtree = apply(out$labelmat.mrtree, 2, FUN=function(x) length(unique(x)))
amri.flat = sapply(1:ncol(out$labelmat.flat), function(i) AMRI(out$labelmat.flat[,i], ref.labels)$amri)
amri.flat = aggregate(amri.flat, by=list(k=ks.flat), FUN=mean)
amri.recon = sapply(1:ncol(out$labelmat.mrtree), function(i) AMRI(out$labelmat.mrtree[,i], ref.labels)$amri)

df = rbind(data.frame(k=amri.flat$k, amri=amri.flat$x, method='Seurat flat'), 
           data.frame(k=ks.mrtree, amri=amri.recon, method='MRtree'))
ggplot2::ggplot(data=df, aes(x=k, y=amri, color=method)) + geom_line() + theme_bw()


# We calcuate the similarity between the initial flat clustering and MRtree clusters across scales. Lower similarity indicates the selected clustering algorithm is not able to generate stabel clusters at the specific resolution. 
stab.out = stability_plot(out)
stab.out$plot

# MRtree with SC3 clustering
ks = 2:10
clusterings.per.k = 1

set.seed(1)
clust.out =  sc_clustering.sc3(exprs=dat$counts, Ks=rep(ks, each=clusterings.per.k), type = 'counts', colData=metadata,
                               biology = F, n_cores = NULL, gene_filter = F, pct_dropout_min=-1)

# Plot the initial clustering reults from SC3

plot_clustree(labelmat=clust.out$sce@colData, prefix='sc3_', suffix = "_clusters",
              ref.labels = ref.labels, plot.ref = F)


#Run MrTree and plot the tree:

out = mrtree(clust.out$sce, prefix='sc3_', suffix = "_clusters", consensus=F, sample.weighted=F)

plot_tree(labelmat=out$labelmat.mrtree, ref.labels=ref.labels, plot.piechart = T,
          node.size = 0.4, tip.label.dist = 10, bottom.margin=30 )


# We evaluate te per-resolution clustering performance via AMRI:

ks.flat = apply(out$labelmat.flat, 2, FUN=function(x) length(unique(x)))
ks.mrtree = apply(out$labelmat.mrtree, 2, FUN=function(x) length(unique(x)))
amri.flat = sapply(1:ncol(out$labelmat.flat), function(i) AMRI(out$labelmat.flat[,i], ref.labels)$amri)
amri.flat = aggregate(amri.flat, by=list(k=ks.flat), FUN=mean)
amri.recon = sapply(1:ncol(out$labelmat.mrtree), function(i) AMRI(out$labelmat.mrtree[,i], ref.labels)$amri)

df = rbind(data.frame(k=amri.flat$k, amri=amri.flat$x, method='SC3 flat'), 
           data.frame(k=ks.mrtree, amri=amri.recon, method='MRtree'))
ggplot2::ggplot(data=df, aes(x=k, y=amri, color=method)) + geom_line() + theme_bw()


# We calcuate the similarity between the initial flat clustering and MRtree clusters across scales.

stab.out = stability_plot(out)
stab.out$plot







message("+-------------------------------------------------------------------------------")
message("+ ONLY ONCE DECIDE ON RESOLUTION: remove clusters smaller than 100 cells        ")
message("+-------------------------------------------------------------------------------")

reso <- "SCT_snn_res.0.8"

# clusters 24 and 25 have 34 and 31 cells respectively.
head(dropseq.integrated@meta.data,2)
clusts_below100 <- as.data.frame(table(dropseq.integrated@meta.data[,reso]))
keep_clusters_over_100 <- clusts_below100[clusts_below100$Freq >= 100,]$Var1

dropseq.integrated <- subset(dropseq.integrated, SCT_snn_res.0.8 %in% keep_clusters_over_100)


message("+-------------------------------------------------------------------------------")
message("+        dot plot of TSC markers                                 ")
message("+-------------------------------------------------------------------------------")

head(dropseq.integrated@meta.data,2)
nrow(dropseq.integrated@meta.data) # after removing clusters < 100 cells :  216547
table(dropseq.integrated@meta.data[reso])


dropseq.integrated@meta.data$Cluster <- as.character(dropseq.integrated@meta.data[,reso])
table(dropseq.integrated@meta.data$Cluster)
table(dropseq.integrated@meta.data$Experiment)

# order by tempora (300g) res1:
#dropseq.integrated@meta.data$Cluster <- factor(dropseq.integrated@meta.data$Cluster , levels = c("12", "6", "18", "1", "19", "10", "23", "16", "2", "14", "11", "5", "15", "0", "9", "17", "8", "13", "4", "20", "3", "21", "7","22"))
# order by tempora (400g) res1:
#dropseq.integrated@meta.data$Cluster <- factor(dropseq.integrated@meta.data$Cluster , levels = c("19", "12", "3", "16", "18", "2", "5", "13", "1", "21", "15", "14", "9", "7", "10", "8", "17", "4", "11", "20", "6", "0"))
# order by tempora (500g) res1:
#dropseq.integrated@meta.data$Cluster <- factor(dropseq.integrated@meta.data$Cluster , levels = c("13", "15", "0", "7", "1", "18", "12", "3", "22", "14","20","21", "5", "2", "11", "6", "9", "8", "16", "17", "19", "10", "4")) # 14...20!
# order by tempora (500g) res0.8:
dropseq.integrated@meta.data$Cluster <- factor(dropseq.integrated@meta.data$Cluster , levels = c("12", "15", "0", "3", "17", "11", "13", "2", "14","20","18", "7", "8", "5", "9", "10", "16", "1", "19", "6", "4")) # 14...20!
dropseq.integrated@meta.data$Cluster <- factor(dropseq.integrated@meta.data$Cluster , levels = c("12", "15", "0", "3", "17", "11", "13", "2", "14","18", "7", "8", "5", "9", "10", "16", "1", "19", "6", "4")) # 14...19!



Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Cluster

min(dropseq.integrated@meta.data$nFeature_RNA)
sum(dropseq.integrated@meta.data$nFeature_RNA > 500)


#DefaultAssay(dropseq.integrated)
DotPlot(dropseq.integrated, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE)

DotPlot(dropseq.integrated, features = c("Tfeb", "Ly6e","Ovol2","Gcm1","Dlx3","Cebpa", "Syna") , cols = c("green", "blue"), assay = "RNA", dot.scale = 10, scale = TRUE)
DotPlot(dropseq.integrated500, features = c("Tfeb", "Ly6e","Ovol2","Gcm1","Dlx3","Cebpa", "Syna") , cols = c("green", "blue"), assay = "RNA", dot.scale = 10, scale = TRUE)


pdf(paste("DotPlot_", Project, "_", "umap.Cluster_no_legend_","res_", reso,".pdf", sep = "_"), width=9, height=5)
par(bg=NA)
DotPlot(dropseq.integrated, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

DefaultAssay(dropseq.integrated) <- "SCT"
FeaturePlot(dropseq.integrated, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1"), ncol = 3)
FeaturePlot(dropseq.integrated, features = c("Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi", "Tfeb", "Ly6e"), ncol = 3)
FeaturePlot(dropseq.integrated, features = c("Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt"), ncol = 3)


message("+-------------------------------------------------------------------------------")
message("+                        Plotting UMAP - custom                                 ")
message("+-------------------------------------------------------------------------------")

sampleTable <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/sampleTable_MERGED_53.csv")

dataset_name <- "Harmony__orig_500g" 

head(dropseq.integrated@meta.data,2)

matrix.umap            <- as.data.frame(Embeddings(object=dropseq.integrated, reduction="umap"))
matrix.umap$Sample_2   <- rownames(matrix.umap) 
matrix.umap$Sample_2   <- gsub("_[AGCTN]{12}$","",           matrix.umap$Sample_2)
sampleTable$Sample_2   <- paste(sampleTable$sampleLabels, sampleTable$sampleLabels, sep = "_")
rownames(sampleTable)  <- sampleTable$Sample_2 
matrix.umap$Sample     <- sampleTable[matrix.umap$Sample_2, ]$sampleLabels
matrix.umap$Experiment <- sampleTable[matrix.umap$Sample_2, ]$CellType
matrix.umap$Batch      <- sampleTable[matrix.umap$Sample_2, ]$Merged_batch
matrix.umap$Project    <- sampleTable[matrix.umap$Sample_2, ]$Project
head(matrix.umap)

cell_cycle             <- dropseq.integrated@meta.data %>% dplyr::select((Phase))  # 
colnames(cell_cycle)   <- c("Phase")
matrix.umap$Phase      <- cell_cycle$Phase

reso                <- "SCT_snn_res.0.8"
clust               <- dropseq.integrated@meta.data %>% dplyr::select((reso))  # 

colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 

head(matrix.umap,2)
unique(matrix.umap$Batch)
unique(matrix.umap$Cluster)


table(matrix.umap$Experiment)
table(matrix.umap$Cluster)

matrix.umap$UMAP_1 <- matrix.umap$UMAP_1 * -1


message("----------------------Cluster annotation -----------------------")

Clusters <- as.numeric(levels(matrix.umap$Cluster))

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


centers <- data.frame(Clusters=Clusters,coord_x=unlist(coord_x), coord_y=unlist(coord_y))
centers <- centers[order(centers$Clusters),]

dat <- matrix.umap

message("----------------------go to script ASSIGN_COLOURS_UMAP.R !!! -----------------------")

gc()

message("----------------------umap.Cluster-----------------------")

table(matrix.umap$Cluster)
matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22")) 



Col_cluster <- c(best.colordf$color)
names(Col_cluster) <- levels(matrix.umap$Cluster)

# for harmony orig 500g (all clusters, nothing removed):::
Col_cluster <- c("#d9f0a3" ,  "#081d58" , "#bfd3e6" , "#fa9fb5" ,  "#88419d", 
                 "#0868ac" ,  "#f768a1" ,  "#4d004b" ,  "#00B3FFFF", "#FF00E6FF",
                 "#FFE500FF", "#004529"  , "#78c679" ,  "#FF0000FF", "#810f7c" ,
                 "#FFFF99" ,  "#00FFB2FF", "#fde0dd" , "#FF9900FF" ,   "#7F00FFFF", 
                 "#d4b9da" ,  "#238443" ,  "#8c96c6" )# , "#4eb3d3"
names(Col_cluster) <- levels(matrix.umap$Cluster)


# for harmony orig 500g reso 0.8 (all clusters, nothing removed):::
Col_cluster <- c( "#FFFF99" , "#d9f0a3" , "#bfd3e6" , "#fa9fb5" ,  "#88419d", 
                  "#0868ac" ,  "#f768a1" ,  "#4d004b" ,  "#00B3FFFF", "#FF00E6FF",
                  "#FFE500FF", "#004529"  , "#78c679" ,  "#FF0000FF", "#810f7c" ,
                  "#081d58" ,  "#00FFB2FF", "#fde0dd" , "#FF9900FF" ,   "#7F00FFFF", 
                  "#d4b9da"  )# , "#4eb3d3"
names(Col_cluster) <- levels(matrix.umap$Cluster)


# for harmony orig 400g (all clusters, nothing removed):::
Col_cluster <- c("#d9f0a3"   "#081d58"   "#FF9900FF" "#fa9fb5"   "#88419d" 
                 "#0868ac"   "#f768a1"   "#4d004b"   "#00B3FFFF" "#FF00E6FF"
                 "#FFE500FF" "#004529"   "#78c679"   "#FF0000FF" "#810f7c" 
                 "#FFFF99"   "#00FFB2FF" "#fde0dd"  "#bfd3e6"   "#7F00FFFF" 
                 "#d4b9da"   "#238443"   "#8c96c6"   "#4eb3d3")
names(Col_cluster) <- levels(matrix.umap$Cluster)




# for harmony orig after removing lcust 24 and 25:::
Col_cluster <- c("#cc99ff", "#0868ac" ,"#8c96c6" ,  "#FF9900FF", "#810f7c",
                 "#004529",   "#4d004b", "#333399"  , "#00B3FFFF" ,"#BF5B17", 
                 "#FF00E6FF" ,"#7F00FFFF", "#FFFF99" ,  "#f768a1"  , "#00FFB2FF",
                 "#78c679"   ,  "#d9f0a3" ,  "#bfd3e6" , "#88419d" ,  "#666666"  , 
                 "#d4b9da"  , "#fa9fb5"  , "#FF0000FF" ,"#fde0dd")
names(Col_cluster) <- levels(matrix.umap$Cluster)



# for harmony orig:::
Col_cluster<-c( "#FF00E6FF", "#FF9900FF", "#FF0000FF", "#238443" , "#fa9fb5" , 
                "#00B3FFFF",   "#88419d" ,  "#810f7c" ,  "#081d58" , "#FFE500FF",
                "#FFFF99"  ,  "#f768a1" ,  "#004529" ,  "#8c96c6" ,  "#4d004b" ,
                "#fde0dd" ,  "#7F00FFFF", "#d9f0a3" ,  "#0868ac"  ,"#d4b9da"  ,
                "#4eb3d3" ,  "#bfd3e6",   "#78c679" ,   "#00FFB2FF")
names(Col_cluster) <- levels(matrix.umap$Cluster)

# cluster anno 500g res 0.8:::
# annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)


# cluster_anno 500g::
#annotate('text', label = paste0('0'), x =12.1268456718008, y = 2.53892382918238, color = 'black', size=8) + annotate('text', label = paste0('1'), x =8.96589071956078, y = 2.81247124014734, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-11.2353831985911, y = 7.46065577803492, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-0.495667541394131, y = 0.39628458706974, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-2.87533049377997, y = 5.54345235167383, color = 'black', size=8) + annotate('text', label = paste0('13'), x =12.4304774543325, y = -4.01739683808447, color = 'black', size=8) + annotate('text', label = paste0('14'), x =-2.62157170567115, y = -3.37125101746679, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.8286569854299, y = 1.33972069559931, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-14.6941301086863, y = 1.7588748627746, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-9.29354016575416, y = 10.9873238497817, color = 'black', size=8) + annotate('text', label = paste0('18'), x =8.7183855315725, y = 11.4531966143691, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.32902710232337, y = 7.34505924521326, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-2.47219722065528, y = -7.59413233460547, color = 'black', size=8) + annotate('text', label = paste0('20'), x =14.9387424727956, y = 4.27432808218836, color = 'black', size=8) + annotate('text', label = paste0('21'), x =-11.7462093094309, y = 1.78895506201624, color = 'black', size=8) + annotate('text', label = paste0('22'), x =-15.5476567009409, y = -6.63459840478064, color = 'black', size=8) + annotate('text', label = paste0('3'), x =5.13639242854516, y = -3.83688679398657, color = 'black', size=8) + annotate('text', label = paste0('4'), x =-11.7088958481272, y = -1.87138107480169, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-7.64920942577918, y = 2.83404561816095, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-10.4104549148997, y = 3.27178034125208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =15.6641557952444, y = 0.442290994658461, color = 'black', size=8) + annotate('text', label = paste0('8'), x =-9.64843766483863, y = -9.01665678681494, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-6.97128121647437, y = -3.51995101632238, color = 'black', size=8)

#    Cluster_Anno 400g:::
#  annotate('text', label = paste0('0'), x =11.2998769033778, y = 1.18818370929158, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-5.03420999161125, y = 3.54964642157948, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.22509095557808, y = 4.2656422101537, color = 'black', size=8) + annotate('text', label = paste0('11'), x =11.2994253385889, y = -2.92744703659618, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3370449792516, y = 3.64214782348073, color = 'black', size=8) + annotate('text', label = paste0('13'), x =1.25368699916481, y = -6.63141126999461, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.58857796081184, y = 3.60154656997121, color = 'black', size=8) + annotate('text', label = paste0('15'), x =3.29137740024208, y = -3.332692817541, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-11.7763779412878, y = -4.43145866760814, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0715372312891, y = -11.8051607168635, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-9.74240806213738, y = -10.522739604803, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-10.0727064859045, y = 0.000445693043829964, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-8.63164451233269, y = -2.89844818482005, color = 'black', size=8) + annotate('text', label = paste0('20'), x =12.1175296057093, y = -5.09075804123485, color = 'black', size=8) + annotate('text', label = paste0('21'), x =1.12023100741981, y = -2.02477283844554, color = 'black', size=8) + annotate('text', label = paste0('22'), x =11.8988244284022, y = -1.93871064552867, color = 'black', size=8) + annotate('text', label = paste0('23'), x =15.1477029074061, y = 8.47610931029713, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-11.6272652398718, y = -2.44086404213512, color = 'black', size=8) + annotate('text', label = paste0('4'), x =9.15782949813484, y = 8.87147168746388, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-14.7408578645361, y = -0.0852867818793039, color = 'black', size=8) + annotate('text', label = paste0('6'), x =8.95425579436897, y = -7.51002498039806, color = 'black', size=8) + annotate('text', label = paste0('7'), x =5.83457825072883, y = -2.81944318184459, color = 'black', size=8) + annotate('text', label = paste0('8'), x =1.80053422340034, y = -0.635730520578263, color = 'black', size=8) + annotate('text', label = paste0('9'), x =2.46752414592384, y = 8.0639865838567, color = 'black', size=8)


#    Cluster_Anno ORIG:::
#annotate('text', label = paste0('0'), x =5.35434635686573, y = 0.898358637808718, color = 'black', size=8) + annotate('text', label = paste0('1'), x =-13, y = 0.138014810501493, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-4.22145811510388, y = -4.04271966743573, color = 'black', size=8) + annotate('text', label = paste0('3'), x =11.5170646529168, y = 1.29690390300647, color = 'black', size=8) + annotate('text', label = paste0('4'), x =8.01382716703113, y = -8.15843135643109, color = 'black', size=8) + annotate('text', label = paste0('5'), x =-0.888398148539696, y = -7.5587819271098, color = 'black', size=8) + annotate('text', label = paste0('6'), x =-11.2329634804756, y = -0.271788959981047, color = 'black', size=8) + annotate('text', label = paste0('7'), x =8.01778276967701, y = -1.46279138851269, color = 'black', size=8) + annotate('text', label = paste0('8'), x =5.11509498166737, y = 6.18444650840656, color = 'black', size=8) + annotate('text', label = paste0('9'), x =-0.816725530150567, y = 0.944908792494692, color = 'black', size=8) + annotate('text', label = paste0('10'), x =-13.4801973481208, y = 4.10048144531146, color = 'black', size=8) + annotate('text', label = paste0('11'), x =1.34165807771381, y = 9.54760044288532, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.0503110070259, y = -5.04408056068524, color = 'black', size=8) + annotate('text', label = paste0('13'), x =14.0734244208306, y = 3.84711270999805, color = 'black', size=8) + annotate('text', label = paste0('14'), x =0.5, y = 2.7, color = 'black', size=8) + annotate('text', label = paste0('15'), x =10.1915961127251, y = 8.15711563300983, color = 'black', size=8) + annotate('text', label = paste0('16'), x =-7.17172137689892, y = 4.2078720874776, color = 'black', size=8) + annotate('text', label = paste0('17'), x =10.0114589552849, y = 14.3196250743856, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-15.5, y = 1.27459519576923, color = 'black', size=8) + annotate('text', label = paste0('19'), x =-7.77777186823193, y = 14.465984041213, color = 'black', size=8) + annotate('text', label = paste0('20'), x =8.7315245490044, y = 3.5356673545827, color = 'black', size=8) + annotate('text', label = paste0('21'), x =8.48272260236439, y = 7.38754838180438, color = 'black', size=8) + annotate('text', label = paste0('22'), x =15.5939398627251, y = -0.180076842785917, color = 'black', size=8) + annotate('text', label = paste0('23'), x =-7.49426738214794, y = 1.94848107767002, color = 'black', size=8) #+ theme(legend.position = "none")




umap.Cluster    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Cluster)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  #ggtitle(paste0(dataset_name, " QC :::  ", reso, " Cluster" )) +
  coord_fixed() + xlab("") + ylab("") +
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=1, size=6))) + 
  theme_classic() + 
  scale_colour_manual("Cluster", values = Col_cluster) +
  #theme(legend.position = "none")  +
  annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)

umap.Cluster
#umap.Cluster_w_legend <- umap.Cluster
#umap.Cluster_no_legend <- umap.Cluster


png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Cluster_w_legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Cluster_no_legend
dev.off()

#saveRDS(umap.Cluster_w_legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.Cluster_w_legend.Rds")
#saveRDS(umap.Cluster_no_legend, "CTR_dscj1_HARMONY_orig_final_FIG_umap.Cluster_no_legend.Rds")

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Cluster_w_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Cluster_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Cluster_no_legend
dev.off()




message("----------------------umap.Batch-----------------------")

#sampleTable_orig <- read.csv("sampleTable_183samples_all_batches_alignrate.csv")
#tbl_for_merging <- read.csv("tbl_for_merging_removedUnsure.csv")


matrix.umap$Batch_orig <- sampleTable_orig[match(matrix.umap$Sample, sampleTable_orig$sampleLabels),]$Batch

umap.batch        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Batch)) + geom_point(alpha=0.4, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  #ggtitle(paste0( dataset_name ,"QC ::: ", reso, " Cell Cycle" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("0B0"="green", "00C"="red",  "0BC"="blue", "A00" = "yellow", "A0C" = "orange", "AB0" = "pink", "ABC" = "black")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) + coord_fixed() + xlab("") + ylab("")+
  theme_classic() #+ theme(legend.position = "none") 
umap.batch

      

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.batch_new_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.batch
dev.off()



message("----------------------umap.CellCycle-----------------------")

umap.cc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Phase)) + geom_point(alpha=0.4, size=0.1) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  #ggtitle(paste0( dataset_name ,"QC ::: ", reso, " Cell Cycle" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("G2M"="red", "G1"="green", "S"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) + coord_fixed() + xlab("") + ylab("")+
  theme_classic() #+ theme(legend.position = "none") 

#umap.cc_w_Legend <- umap.cc
#umap.cc_no_Legend <- umap.cc


png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.cc_w_Legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.cc_no_Legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.cc_w_Legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.cc_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.cc_no_Legend
dev.off()


message("----------------------umap.Experiment-----------------------")

umap.Experiment    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.2, size=0.3) +
  #ggtitle(paste0("Resolution ", reso, " :: Experiment " )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Experiment, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="black", "I.1"="yellow", "R.1"="yellow2", 
                                     "I.4"="orangered", "R.4"="red", "I.24"="lightblue2", "R.24"="lightblue4",
                                     "I.36"="blue", "R.36"="darkblue","I.48"="purple", "R.48"="orchid"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  theme_classic() + coord_fixed() +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(title=element_text(size=6), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10))# + theme(legend.position = "none") 

#umap.Expt_w_Legend <- umap.Experiment
#umap.Expt_no_legend <- umap.Experiment


png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_w_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Expt_w_legend
dev.off()

png(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_no_legend_","res_", reso,".png"), width=1000, height=900, type="cairo")
par(bg=NA)
umap.Expt_no_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_w_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Expt_w_legend
dev.off()

pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "umap.Expt_no_legend_","res_", reso,".pdf"), width=10, height=9)
par(bg=NA)
umap.Expt_no_legend
dev.off()




Cluster_Anno

message("----------------------umap.t.0-----------------------")

matrix.umap$Mixer <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_timepoints")
table(matrix.umap$Mixer)

umap.t0    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Mixer)) +
  geom_point(alpha=0.3, size=0.3) +
  #ggtitle(paste0("Resolution ", reso, " T0 in Batches" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("other_timepoints"="lightgrey", "0.0" = "black")) +
  coord_fixed() +
  xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  #guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  theme(legend.position = "none") +
  annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)






message("----------------------umap.t.1-----------------------")

umap.Experiment1    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="blue", "R.1"="red", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)






message("----------------------umap.t.4-----------------------")

umap.Experiment4    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="blue", "R.4"="red", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)






message("----------------------umap.t.24-----------------------")

umap.Experiment24    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="blue", "R.24"="red",
                                     "I.36"="grey", "R.36"="grey","I.48"="grey", "R.48"="grey"), 
                      limits=c('0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))  + 
  annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)





message("----------------------umap.t.36-----------------------")

umap.Experiment36    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="blue", "R.36"="red","I.48"="grey", "R.48"="grey"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)



message("----------------------umap.t.48-----------------------")

umap.Experiment48    <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=Experiment)) +
  geom_point(alpha=0.1, size=0.1) +
  #ggtitle(paste0("Resolution ", reso, " Experiment" )) +
  geom_density_2d(data=matrix.umap, aes(x=UMAP_1, y=UMAP_2, group=Cluster, alpha=..level..), 
                  linetype='dashed', colour="black", bins=3, alpha=0.5) +
  scale_colour_manual("", values = c("0.0"="grey", "I.1"="grey", "R.1"="grey", 
                                     "I.4"="grey", "R.4"="grey", "I.24"="grey", "R.24"="grey",
                                     "I.36"="grey", "R.36"="grey","I.48"="blue", "R.48"="red"), 
                      limits=c( '0.0', 'I.1', "R.1", 'I.4', "R.4", 'I.24', "R.24", 'I.36', "R.36", 'I.48', "R.48")) + 
  coord_fixed() + xlab("") + ylab("") + theme_classic() + 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  annotate('text', label = paste0('0'), x =-11.7201665183584, y = 2.53170069990992, color = 'black', size=8) + annotate('text', label = paste0('1'), x =9.12361733707984, y = -8.75345507325293, color = 'black', size=8) + annotate('text', label = paste0('2'), x =-5.12661011424462, y = -3.74864235581518, color = 'black', size=8) + annotate('text', label = paste0('3'), x =-8.96201784816186, y = 2.71836289702295, color = 'black', size=8) + annotate('text', label = paste0('4'), x =11.3610984543284, y = -0.753229886904964, color = 'black', size=8) + annotate('text', label = paste0('5'), x =10.91468159947, y = 2.97537264166712, color = 'black', size=8) + annotate('text', label = paste0('6'), x =11.1705972412546, y = 7.42909583388208, color = 'black', size=8) + annotate('text', label = paste0('7'), x =0.864693866105871, y = -5.21964111985327, color = 'black', size=8) + annotate('text', label = paste0('8'), x =7.72890727314551, y = 2.89469608603357, color = 'black', size=8) + annotate('text', label = paste0('9'), x =1.57538645062049, y = 0.517522923036804, color = 'black', size=8) + annotate('text', label = paste0('10'), x =6.63360445293983, y = -3.81958547295691, color = 'black', size=8) + annotate('text', label = paste0('11'), x =-15.7878101607839, y = 0.249652101978054, color = 'black', size=8) + annotate('text', label = paste0('12'), x =-12.3741081496755, y = -4.09128561677099, color = 'black', size=8) + annotate('text', label = paste0('13'), x =2.90870849880775, y = 5.70920333205103, color = 'black', size=8) + annotate('text', label = paste0('14'), x =2.61351244244178, y = -3.365908414737, color = 'black', size=8) + annotate('text', label = paste0('15'), x =-10.7590120574514, y = 2.37340518771051, color = 'black', size=8) + annotate('text', label = paste0('16'), x =9.29303948673804, y = 10.9854613238418, color = 'black', size=8) + annotate('text', label = paste0('17'), x =-8.71827013698022, y = 11.4546938830459, color = 'black', size=8) + annotate('text', label = paste0('18'), x =-14.9276278754751, y = 4.2959899598205, color = 'black', size=8) + annotate('text', label = paste0('19'), x =7.33367841038306, y = 7.32581648169397, color = 'black', size=8) + annotate('text', label = paste0('20'), x =15.5476567009409, y = -6.63459840478064, color = 'black', size=8)




png(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "umap__Expt_grid_eachTimePoint", reso, ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
plot_grid(umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 2, align = "hv")
dev.off()


#####
list_umaps_for_plotting_noLegend <- list(umap.t0=umap.t0, 
                                         umap.Experiment1=umap.Experiment1, 
                                         umap.Experiment4=umap.Experiment4, 
                                         umap.Experiment24=umap.Experiment24, 
                                         umap.Experiment36=umap.Experiment36, 
                                         umap.Experiment48=umap.Experiment48)
#saveRDS(list_umaps_for_plotting_noLegend,"list_umaps_for_plotting_noLegend_500g_reso0.8.Rds")
list_umaps_for_plotting_2 <- list(umap.cc_no_Legend=umap.cc_no_Legend, 
                                  umap.cc_w_Legend =umap.cc_w_Legend, 
                                  umap.Cluster_no_legend  =umap.Cluster_no_legend, 
                                  umap.Cluster_w_legend=umap.Cluster_w_legend, 
                                  umap.Expt_no_legend =umap.Expt_no_legend, 
                                  umap.Expt_w_Legend=umap.Expt_w_Legend )
#saveRDS(list_umaps_for_plotting_2,"list_umaps_for_plotting2_cc_clust_expt_500g_reso0.8.Rds")


#list_umaps_for_plotting_2        <- readRDS("list_umaps_for_plotting2_cc_clust_expt.Rds")
##list_umaps_for_plotting_noLegend <- readRDS("list_umaps_for_plotting_noLegend.Rds")
plt_pseudotime <- readRDS("MONOCLE_SCT/harmony-orig_monocle2_plt_pseudotime_loopT.Rds")

names(list_umaps_for_plotting_noLegend)
names(list_umaps_for_plotting_2)

pdf(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "9x_GRID_PAPER_FIGURE", reso, ".pdf", sep="_"), width=20,height=18)
#png(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "9x_GRID_PAPER_FIGURE", reso, ".png", sep="_"), width=2500,height=2000, type = "cairo")
par(bg=NA)
plot_grid( list_umaps_for_plotting_2[["umap.Cluster_no_legend"]], 
           list_umaps_for_plotting_2[["umap.cc_no_Legend"]] + theme(legend.position = "none")+ xlab("") + ylab("") + coord_fixed() ,
           NA , 
           list_umaps_for_plotting_noLegend[["umap.t0"]], 
           list_umaps_for_plotting_noLegend[["umap.Experiment1"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment4"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment24"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment36"]],
           list_umaps_for_plotting_noLegend[["umap.Experiment48"]], ncol = 3, align = "vh")
dev.off()
# plt_pseudotime + theme(legend.position = "none") + xlab("") + ylab("")+ coord_fixed() , 


gc()

#need to calculate pseudotime!!! 
png(paste( "UMAP_MergedRealigned_51_", Project, "ggplot", "GRID_PAPER_FIGURE", reso, ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
plot_grid(umap.Cluster_no_legend, umap.cc_no_Legend, umap.pseudotime, list_umaps_for_plotting_noLegend["umap.t0"], list_umaps_for_plotting_noLegend["umap.Experiment1"],list_umaps_for_plotting_noLegend["umap.Experiment4"],list_umaps_for_plotting_noLegend["umap.Experiment24"],list_umaps_for_plotting_noLegend["umap.Experiment36"],list_umaps_for_plotting_noLegend["umap.Experiment48"], ncol = 3, nrow = 3, align = "hv")
dev.off()




pdf(paste0("UMAP_MergedRealigned51_", Project,  "_ggplot_", "GRID_PAPER_FIGURE","res_", reso,".pdf"), width=20, height=10)
par(bg=NA)
plot_grid(umap.Cluster_no_legend, umap.cc_no_Legend, umap.pseudotime, umap.t0, umap.Experiment1, umap.Experiment4, umap.Experiment24, umap.Experiment36, umap.Experiment48, ncol = 3, nrow = 3, align = "hv")
dev.off()




message("----------------------umap.QC_1_sample-----------------------")

table(dropseq.integrated@meta.data[,c(5,8)])

matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0", "0.0", "other_samples")
matrix.umap$BatchSample <- paste( matrix.umap$Batch, matrix.umap$Sample, sep = "_")
unique(matrix.umap[matrix.umap$QC == "0.0",]$BatchSample)
# "0B0_Batch0011_15988_N704" "0BC_Dups_setA1_reps_N724" "A00_Dups_setA2_reps_N719" "A00_Dups_setB_reps_N704"  "00C_Batch0004_9617_N722" 
# "00C_Batch0004_9617_N726"  "00C_Batch0004_9617_N727"  "00C_Batch0004_9617_N720"  "00C_Batch0004_9617_N719"  "00C_Batch0004_9617_N718" 
#  "00C_Batch0004_9617_N729" 

matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "0B0_Batch0011_15988_N704", "0B0_Batch0011_15988_N704", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "0BC_Dups_setA1_reps_N724", "0BC_Dups_setA1_reps_N724", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "A00_Dups_setA2_reps_N719", "A00_Dups_setA2_reps_N719", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "A00_Dups_setB_reps_N704", "A00_Dups_setB_reps_N704", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N722", "00C_Batch0004_9617_N722", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N726", "00C_Batch0004_9617_N726", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N727", "00C_Batch0004_9617_N727", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N720", "00C_Batch0004_9617_N720", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N719", "00C_Batch0004_9617_N719", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N718", "00C_Batch0004_9617_N718", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "0.0" & matrix.umap$BatchSample == "00C_Batch0004_9617_N729", "00C_Batch0004_9617_N729", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.3, size=0.1) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.0 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  #scale_colour_manual("", values = c("other_samples"="lightgrey", "A00_t0" = "blue1", "00C_t0"="red", "0B0_t0" = "orange", "0BC_t0" = "green")) +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "0B0_Batch0011_15988_N704" = "blue1", "0BC_Dups_setA1_reps_N724"="red", "A00_Dups_setA2_reps_N719" = "orange", "A00_Dups_setB_reps_N704" = "green", "00C_Batch0004_9617_N722"= "khaki3", "00C_Batch0004_9617_N726" = "brown", "00C_Batch0004_9617_N727" = "navyblue", "00C_Batch0004_9617_N720"= "pink", "00C_Batch0004_9617_N719"= "violet", "00C_Batch0004_9617_N718"= "turquoise")) + # , "00C_Batch0004_9617_N729" = "black"
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() #+ xlim(c(6,7)) + ylim(c(-3, -1.8)) # + xlim(c(3,9)) + ylim(c(-3.5, 2.5))
umap.qc


matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1", "I.1", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "0B0", "0B0_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "A00", "A00_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "A0C", "A0C_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "ABC", "ABC_I.1", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.1" & matrix.umap$Batch == "AB0", "AB0_I.1", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.1 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "A00_I.1" = "blue1", "A0C_I.1"="red", "0B0_I.1" = "orange", "ABC_I.1" = "green", "AB0_I.1"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24", "I.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$Batch == "ABC", "ABC_I.24", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$Batch == "AB0", "AB0_I.24", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC T.24 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_I.24" = "green", "AB0_I.24"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc



matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4", "R.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "ABC", "ABC_R.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "AB0", "AB0_R.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$Batch == "A0C", "A0C_R.4", matrix.umap$QC)
unique(matrix.umap$QC)

umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC r.4 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_R.4" = "green", "AB0_R.4"="purple", "A0C_R.4"="blue")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4", "I.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "ABC", "ABC_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "AB0", "AB0_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "A0C", "A0C_I.4", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.4" & matrix.umap$Batch == "A00", "A00_I.4", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC i.4 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_I.4" = "green", "AB0_I.4"="purple", "A0C_I.4"="blue", "A00_I.4" = "orange")) +   guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6))) + 
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=14), text=element_text(size=14),  axis.text.x=element_text(size=14), axis.text.y=element_text(size=14)) +
  theme_classic() 
umap.qc




matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24", "R.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24" & matrix.umap$Batch == "ABC", "ABC_R.24", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.24" & matrix.umap$Batch == "A00", "A00_R.24", matrix.umap$QC)
#matrix.umap$QC <- ifelse(matrix.umap$Sample_2 == "Batch0004_9617_N721", "Batch0004_9617_N721", matrix.umap$QC)
unique(matrix.umap$QC)


umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC r.24 batches" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_R.24" = "green", "A00_R.24"="purple")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc


unique(matrix.umap[matrix.umap$Experiment == "R.36",]$BatchSample)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36", "R.36", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R4_reps_N712", "0B0_Dups_36R4_reps_N712", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R2_reps_N707", "0B0_Dups_36R2_reps_N707", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R1_reps_N706", "0B0_Dups_36R1_reps_N706", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.36" & matrix.umap$BatchSample == "0B0_Dups_36R3_reps_N710", "0B0_Dups_36R3_reps_N710", matrix.umap$QC)

umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC R.36 samples" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "0B0_Dups_36R4_reps_N712" = "green", "0B0_Dups_36R2_reps_N707"="purple", "0B0_Dups_36R1_reps_N706"="blue", "0B0_Dups_36R3_reps_N710"="red")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc



unique(matrix.umap[matrix.umap$Experiment == "R.4",]$BatchSample)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4", "R.4", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "ABC_Dups_5A_reps_N726", "ABC_Dups_5A_reps_N726", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "AB0_Dups_6C_reps_N720", "AB0_Dups_6C_reps_N720", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "ABC_Dups_6B_reps_N715", "ABC_Dups_6B_reps_N715", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "R.4" & matrix.umap$BatchSample == "A0C_Dups_6D_reps_N720", "A0C_Dups_6D_reps_N720", matrix.umap$QC)

umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC R.4 samples" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_Dups_5A_reps_N726" = "green", "AB0_Dups_6C_reps_N720"="purple", "ABC_Dups_6B_reps_N715"="blue", "A0C_Dups_6D_reps_N720"="red")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc



unique(matrix.umap[matrix.umap$Experiment == "I.24",]$BatchSample)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24", "I.24", "other_samples")
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "ABC_Dups_8D_reps_N722", "ABC_Dups_8D_reps_N722", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "ABC_Dups_7A_reps_N728", "ABC_Dups_7A_reps_N728", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "AB0_Dups_8C_reps_N706", "AB0_Dups_8C_reps_N706", matrix.umap$QC)
matrix.umap$QC <- ifelse(matrix.umap$Experiment == "I.24" & matrix.umap$BatchSample == "ABC_Dups_8B_reps_N710", "ABC_Dups_8B_reps_N710", matrix.umap$QC)

umap.qc        <- ggplot(matrix.umap, aes(x=UMAP_1, y=UMAP_2, colour=QC)) + geom_point(alpha=0.4, size=0.15) +
  ggtitle(paste0( dataset_name ,"QC ::: ", reso, " QC I.24 samples" )) +  coord_fixed() + xlab("") + ylab("") +
  scale_colour_manual("", values = c("other_samples"="lightgrey", "ABC_Dups_8D_reps_N722" = "green", "ABC_Dups_7A_reps_N728"="purple", "AB0_Dups_8C_reps_N706"="blue", "ABC_Dups_8B_reps_N710"="red")) +
  guides(colour = guide_legend(reverse=T, override.aes = list( alpha=0.5, size=6)))+ 
  theme(title=element_text(size=10), text=element_text(size=10), axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
  theme_classic() 
umap.qc










message("+-------------------------------------------------------------------------------")
message("+                 plot previously identified genes                              ")
message("+-------------------------------------------------------------------------------")

# top QC check genes: 
FeaturePlot(dropseq.integrated, features = c("Elf5", "Eomes", "Sox2","Ascl2","Gcm1","Ovol2", "Anxa1" , "Plac1", "Uba6")) # "Anxa1,"Pparg" "Cdx2"

png(paste( "FeaturePlt_UMAP_MergedRealigned_52_int_", Project, "theta4_ref0B0","ggplot",  ".png", sep="_"), width=2000,height=1000, type = "cairo")
par(bg=NA)
FeaturePlot(dropseq.integrated, features = c("Elf5", "Eomes", "Sox2","Ascl2","Gcm1","Ovol2", "Anxa1" , "Plac1", "Uba6")) # "Anxa1,"Pparg" "Cdx2"
dev.off()



FeaturePlot(dropseq.integrated, features = c("Elf5", "Cdx2","Eomes","Tead4","Sox2","Ascl2","Gcm1","Ovol2"))
# LP ("Gcm1","Ovol2"), EPC = Ascl2


FeaturePlot(dropseq.integrated, features = c("Usf2", "Nfya", "Yy1", "Maff", "Gcm1","Ovol2", "Pparg","Ghrl1")) # LP
FeaturePlot(dropseq.integrated, features = c("Vezf1", "Gata1", "Klf7", "Bdp1","Ascl2")) # EPC
FeaturePlot(dropseq.integrated, features = c("Bhlhe40", "Smarca4", "Myc", "Foxo3","Sp3","Srebf2","Hdac6","Tead3", "Ppard")) # TFs
# other TFs : "E2f8","","Grhl1","","Foxo4","Hdac2","","Junb","Nfkb2","Creb3l2","Pou2f1","Mycn","","Max","","","E2f7","Fos","Elf2"


VlnPlot(dropseq.integrated, features = c("Lgals3", "Krt18"), assay = "SCT", log = TRUE)
FeaturePlot(dropseq.integrated, features = c("Flt1", "Anxa1", "Uba6", "Bnip3","Tceb2", "Plac1")) # earlier markrs




# T0 identified genes: 
FeaturePlot(dropseq.integrated, features = c("Krt18", "Lgals3", "Malat1","Rhox6","Atrx","Cenpf", "Ubb" , "Nars", "Ncl")) # "Anxa1,"Pparg" "Cdx2"
FeaturePlot(dropseq.integrated, features = c("Rhox9", "Rpl37", "Krt8","Crip1","Igf2","Cdkn1c", "Npm1" , "Hand1", "Rpl41")) # "Anxa1,"Pparg" "Cdx2"




message("+-------------------------------------------------------------------------------")
message("+                        FIND MARKERS for res 0.8                               ")
message("+-------------------------------------------------------------------------------")

Project <- "CTR_dscj1_HARMONY_500g_res0.8_markers"

markerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/Harmony.orig_500g_Markers"
setwd(markerDir)

# set up variables:
l2fc_cutoff <- 0.25
reso <- "res.0.8"
last_cluster <- length(unique(dropseq.integrated@meta.data[, paste0("SCT_snn_", reso)]))-1

head(dropseq.integrated@meta.data,2)
Idents(dropseq.integrated) <-  dropseq.integrated@meta.data[,paste0("SCT_snn_", reso)]

cluster_numbers <- c(0:(  length(unique(dropseq.integrated@meta.data[, paste0("SCT_snn_", reso)]))-1 ))
marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = length(cluster_numbers), nrow = length(cluster_numbers) ))

marker_tbl <- dplyr::as_tibble(matrix(NA, ncol = last_cluster+1, nrow = last_cluster+1 ))
rownames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20") 
colnames(marker_tbl) <- c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20") 

table(dropseq.integrated@meta.data[, paste0("SCT_snn_", reso)])



dataset_name <- "CTR_dscj1_HARMONY.orig_500g_res0.8_markers"

library(purrr)

findMarkers_for_cluster_pair <- function(i){
  if (j == i){
    print("same cluster")
  } else {
    tmp_markers <- FindMarkers(dropseq.integrated, ident.1 = j, ident.2= i, logfc.threshold = l2fc_cutoff, verbose = TRUE)
    print(length(unique(row.names(tmp_markers))))
    return(tmp_markers)
  }
}
safely_findMarkers = safely( findMarkers_for_cluster_pair )


result_list <- list()
for (n in 0:last_cluster){
  cluster_numbers <- c(n:last_cluster)
  j <- as.character(n)
  result = map(cluster_numbers, safely_findMarkers)
  saveRDS(result, paste0(dataset_name, reso,"_result_", j, "_l2fc" , l2fc_cutoff ,  ".rds"))
  result_list[[n+1]] <- result
}
names(result_list) <- c("result_00", "result_01","result_02","result_03","result_04","result_05","result_06","result_07","result_08","result_09","result_10","result_11","result_12","result_13","result_14","result_15","result_16","result_17","result_18","result_19","result_20") #,"result_21","result_22")
saveRDS(result_list, "result_list_harmony.orig_500g_res0.8.Rds")






result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/Harmony.orig_500g_Markers/result_list_harmony.orig_500g_res0.8.Rds")

marker_tbl  <- readRDS( "marker_tbl.Rds")



message("+-------------------------------------------------------------------------------")
message("+                         MARKERS for res 0.8                                   ")
message("+-------------------------------------------------------------------------------")

markerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/Harmony.orig_500g_Markers"
setwd(markerDir)

reso <- "res.0.8"
l2fc_cutoff <- 0.6
#result_list <- readRDS("HARMONY_together_MARKERS/result_list_harmony_together_res.1.0_l2fc0.6.Rds")
#marker_tbl <- readRDS( "HARMONY_together_MARKERS/marker_tbl.Rds")
#result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/result_list_harmony_together_res.1.0_l2fc0.25.Rds")
#marker_tbl <- readRDS( "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/marker_tbl_res.1.0_l2fc0.25.Rds")
last_cluster <- 20


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
  length(marker_tbl_list[[z]])
}
marker_tbl_list[[1]]


marker_tbl <- as.data.frame(marker_tbl)
length(marker_tbl_list)
for (z in seq_along(marker_tbl_list)){
  marker_tbl[z,] <- c( rep(NA, (last_cluster+1)-length(marker_tbl_list[[z]]) ), marker_tbl_list[[z]])
}
marker_tbl

length(marker_tbl_list[[1]])
length(marker_tbl_list[[2]])
length(marker_tbl_list[[3]])
length(marker_tbl_list[[4]])

max(marker_tbl, na.rm = T) # 78



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
      #tmp_df <- unique(merge(tmp_df, ensEMBL2id, by.x = "row.names", by.y = "external_gene_name", all.x = TRUE))
      tmp_df <- subset(tmp_df, tmp_df$p_val_adj < 0.05)
      tmp_df <- subset(tmp_df, abs(tmp_df$avg_logFC) > l2fc_cutoff)
      write.csv(tmp_df, paste0("CTR_dscj1_HARMONY.orig_500g_" , reso, "_l2fc", l2fc_cutoff,"_Markers_tbl__clusters_", name_y, ".vs.", name_x, ".csv"))
      print(nrow(y[[x]])) }
  }
}

#save_res_tables(x=1, y=result_01)
#save_res_tables(x=4, y=result_01)


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
f3 = colorRamp2( c(0, 3, 6, 9,12), c("blue4","darkorchid4","maroon3","orchid1","white" ), space = "RGB") 

ht3 = Heatmap(as.matrix(marker_tbl),  col = f3, row_title = "", column_title = paste0("Markers between cluster pairs_", reso,  "_absl2fc", l2fc_cutoff), show_row_names = TRUE, heatmap_legend_param = list(title = "Number of markers", legend_height = unit(8, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left",  show_heatmap_legend = TRUE) #  width = unit(10, "cm"),
ht3


pdf(paste("Fig__Pairwise_MARKERS", Project, "ComplexHeatmap",  reso, "_l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=8, height=8)
par(bg=NA)
draw(ht3, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()

#write.csv(marker_tbl, "CTR_dscj1_harmony.orig_500g__MARKERS___Pairwise_cluster_marker_tbl__res.1_l2f0.25.csv", quote = FALSE)
#saveRDS(marker_tbl, "harmony.orig_500g__marker_tbl_res.1.0_l2fc0.25.Rds")




message("--------------------------------------------------------------------------------")
message("+          Now caculate FIND TOP MARKERS for all comparisons                    ")
message("+-------------------------------------------------------------------------------")

library(biomaRt)
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene_id'), mart = ensembl)          



#result_list <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/Harmony.orig_500g_Markers_GO/NewAlign51_SCT_HARMONY_500g_res0.8_GO_res_list_l2fc0.25.Rds")


l2fc_cutoff <- 0.6
marker_genes <- "genes"

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
      print(tmp_df$Row.names) 
      return(tmp_df$Row.names)
    }
  }
}

marker_gene_list <- list()
for (z in seq_along(result_list)){
  gene_list <- list()
  res <- result_list[[z]]
  res_name <- names(result_list)[[z]]
  for (w in 1:length(res)){
    gene_list[[w]] <- save_res_tables(w, res)
  }
  marker_gene_list[[z]] <- gene_list
}


#marker_genes_l2fc1.5 <- unique(unlist(unlist(marker_gene_list)))
#marker_genes_l2fc1 <- unique(unlist(unlist(marker_gene_list)))
marker_genes_l2fc0.6 <- unique(unlist(unlist(marker_gene_list)))
marker_genes_l2fc0.6 <- marker_genes_l2fc0.6[-1]

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Cluster

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[1:40] , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[41:80] , cols = c("green", "blue"), assay = "SCT",  dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[81:120] , cols = c("green", "blue"), assay = "SCT",  dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[121:160] , cols = c("green", "blue"), assay = "SCT",  dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[161:200] , cols = c("green", "blue"), assay = "SCT",  dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[201:240] , cols = c("green", "blue"), assay = "SCT",  dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[241:280] , cols = c("green", "blue"), assay = "SCT",  dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(dropseq.integrated, features = marker_genes_l2fc0.6[281:305] , cols = c("green", "blue"), assay = "SCT",  dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#mat_clusterMeans <- GetAssayData(dropseq.integrated, slot = "scale.data")
#mat_clusterMeans <- mat_clusterMeans[rownames(mat_clusterMeans) %in% marker_genes_l2fc0.6,]
#dim(mat_clusterMeans) # [1]    231 136291
#means_clust <- 
#dim(dropseq.integrated[["SCT"]]@scale.data) # 3000 136291
mat_clusterMeans <- AverageExpression(dropseq.integrated, assays = "SCT", slot = "scale.data")
mat_clusterMeans_SCT <- mat_clusterMeans[[1]][rownames(mat_clusterMeans[[1]]) %in% marker_genes_l2fc0.6,]
mat_clusterMeans_SCT[1:5, 1:5]
mat_clusterMeans_SCT <- mat_clusterMeans_SCT - rowMeans(mat_clusterMeans_SCT)
min_val <- min(mat_clusterMeans_SCT)
max_val <- max(mat_clusterMeans_SCT)

clust_order <- c("12", "15", "0", "3", "17", "11", "13", "2", "14","20","18", "7", "8", "5", "9", "10", "16", "1", "19", "6", "4")
mat_clusterMeans_SCT <- mat_clusterMeans_SCT[, match(clust_order, colnames(mat_clusterMeans_SCT))]

# for slot: data:
f1 = colorRamp2( c(-10, -3, 0, 3, 10), c("#173715", "green3", "grey95", "#bf7acd", "#4b0e81"), space = "LAB") 
# for scale.data:
f1 = colorRamp2( c(min_val, 0, 3, 10), c("green3", "grey95", "#bf7acd", "#4b0e81"), space = "LAB") 
lgd1 = Legend(col_fun = f1, title = "Mean cluster expression", at = c(min_val, max_val )  )

ht1 = Heatmap(t(mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "Marker genes (abs l2fc > 0.6)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE ,  row_names_side ="left",  row_title_rot = 90, column_title_rot = 0 , na_col = "lightgrey", row_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "bottom", column_names_side = "top", row_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1
ht1

pdf(paste(Project, "Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_clusterOrderedTempora.pdf", sep="_"), onefile=FALSE, width=7, height=45) 
par(bg=NA)
draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


ht1_t = Heatmap(t(mat_clusterMeans_SCT[,]),  col = f1, name = "  ",  row_title = " ", column_title = "Marker genes (abs l2fc > 0.6)", show_row_names = TRUE, heatmap_legend_param = list(title = "  ", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = FALSE ,  row_names_side ="left",  row_title_rot = 90, column_title_rot = 0 , na_col = "lightgrey", column_names_gp = gpar(fontface = "italic"), row_dend_side = "right", column_dend_side = "top", column_names_side = "bottom",  row_dend_reorder = TRUE, column_dend_reorder = TRUE) # width = unit(15, "cm"), split = Split$condition, col = col_mat_vst1
ht1_t

pdf(paste(Project, "t_Htmap_ALL_marker_genes_l2fc0.6", "_meanCentr", "scale.data_clusterOrderedTempora.pdf", sep="_"), onefile=FALSE, width=50, height=7) 
par(bg=NA)
draw(ht1_t, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


#library(seriation)
#o1 = seriate(dist(mat_clusterMeans_SCT), method = "TSP")
#o2 = seriate(dist(t(mat_clusterMeans_SCT)), method = "TSP")
#Heatmap(mat, name = "mat", row_order = get_order(o1), column_order = get_order(o2),
#        column_title = "seriation from the distance matrix")




avgexp = AverageExpression(dropseq.integrated, return.seurat = T, add.ident = 'Experiment')
avgclust = AverageExpression(dropseq.integrated, return.seurat = T)

Seurat::DoHeatmap(avgclust, features = marker_genes_l2fc0.6[41:80])

# try... devtools::install_github("Simon-Leonard/FlexDotPlot")


DotPlot(dropseq.integrated, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt") , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



message("--------------------------------------------------------------------------------")
message("+                     GO  using enrichR                               ")
message("+-------------------------------------------------------------------------------")

library("enrichR")
enrichR_DB <- as.data.frame(listEnrichrDbs())

#l2fc_cutoff <- 0.6
l2fc_cutoff <- 0.25

#markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_together_MARKERS/HARMONY_together_MARKERS_GO_res1.0_l2fc0.6"
markerGODir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/Harmony.orig_500g_Markers_GO"
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
        write.csv(enrichR_RES, paste0(Project, "_enrichR_", db, "_l2fc", l2fc_cutoff,"_", comparison_name, ".csv"))
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


FeaturePlot(dropseq.integrated, features = c("Saa3", "Ccl2", "Gata2", "Gata3", "Sox2", "Oct4", "Nanog", "Dab2"))
FeaturePlot(dropseq.integrated, features = c( "Klf17", "Tfap2c", "Slc28a3", "Adap2", "Igfbp3", "Bamb1", "Havcr1"))
# Cga

#### DECIDE ON PARAMETERS:::
#l2fc_cutoff <- 0.6
l2fc_cutoff <- 0.25
database <- "KEGG_2019_Mouse"
# WikiPathways_2019_Mouse, GO_Biological_Process_2018, GO_Cellular_Component_2018, GO_Molecular_Function_2018,  BioCarta_2016, KEGG_2019_Mouse

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

#WP_matrix_l2fc0.25 <- GO_matrix
#GOBP_matrix_l2fc0.25 <- GO_matrix
#GOCC_matrix_l2fc0.25 <- GO_matrix
#GOMF_matrix_l2fc0.25 <- GO_matrix
#BioCarta_matrix_l2fc0.25 <- GO_matrix
#Kegg_matrix_l2fc0.25 <- GO_matrix

#WP_matrix_l2fc0.6  <- GO_matrix
#GOBP_matrix_l2fc0.6  <- GO_matrix
#GOCC_matrix_l2fc0.6  <- GO_matrix
#GOMF_matrix_l2fc0.6  <- GO_matrix
#BioCarta_matrix_l2fc0.6  <- GO_matrix
#Kegg_matrix_l2fc0.6 <- GO_matrix


GO_res_list <- list(GOBP_matrix_l2fc0.25=GOBP_matrix_l2fc0.25,
                    GOCC_matrix_l2fc0.25=GOCC_matrix_l2fc0.25,
                    GOMF_matrix_l2fc0.25=GOMF_matrix_l2fc0.25,
                    BioCarta_matrix_l2fc0.25=BioCarta_matrix_l2fc0.25,
                    WP_matrix_l2fc0.25=WP_matrix_l2fc0.25,
                    Kegg_matrix_l2fc0.25=Kegg_matrix_l2fc0.25)

#GO_res_list <- list(GOBP_matrix_l2fc0.6=GOBP_matrix_l2fc0.6,
#                    GOCC_matrix_l2fc0.6=GOCC_matrix_l2fc0.6,
#                    GOMF_matrix_l2fc0.6=GOMF_matrix_l2fc0.6,
#                    BioCarta_matrix_l2fc0.6=BioCarta_matrix_l2fc0.6,
#                    WP_matrix_l2fc0.6=WP_matrix_l2fc0.6,
#                    Kegg_matrix_l2fc0.6=Kegg_matrix_l2fc0.6)

#saveRDS(GO_res_list, "NewAlign51_SCT_HARMONY_500g_res0.8_GO_res_list_l2fc0.25.Rds")
GO_res_list <- readRDS("NewAlign51_SCT_HARMONY_500g_res0.8_GO_res_list_l2fc0.25.Rds")

names(GO_res_list)

for(i in seq_along(GO_res_list)){
  GO_res_list[[i]] <- GO_res_list[[i]][-c(1,2),]
}
head(GO_res_list[[1]],2)

for(i in seq_along(GO_res_list)){
  write.csv(GO_res_list[[i]], paste( "NewAlign51_SCT_HARMONY_500g__res_0.8_", names(GO_res_list)[[i]], ".csv", sep = "_"),  quote = FALSE)
}




#names(GO_res_list)[6] <-"Kegg_matrix_l2fc0.6" 
#names(GO_res_list)[12] <-"Kegg_matrix_l2fc0.25"

GO_matrix_for_plotting <- GO_res_list[[1]]
GO_res_name <- names(GO_res_list)[[1]]

for(i in seq_along(GO_res_list)){
  #for(i in c(7,8)){
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
  rownames(GO_matrix2) <- gsub( " Homo sapiens h.*" ,  "" , rownames(GO_matrix2))
  
  dim(GO_matrix2)
  f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
  ht1 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , row_names_side ="left", width = unit(ncol(GO_matrix2), "cm")) # width = unit(140, "cm"),
  ht1
  
  pdf(paste( Project, "ComplexHeatmap", "Fig__ALL_Pairwise_", GO_res_name, reso, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix2)/2+5), height=nrow(GO_matrix2)/2)
  par(bg=NA)
  draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
}
library(seriation)






message("--------------------------------------------------------------------------------")
message("+                 which cluster is T0?????                              ")
message("+-------------------------------------------------------------------------------")

head(dropseq.integrated@meta.data,2)
table(dropseq.integrated@meta.data[,c("SCT_snn_res.0.8", "Experiment")])


T0_clust <- 12 # together

GO_matrix_for_plotting <- GO_res_list[[4]]
GO_res_name <- names(GO_res_list)[[4]]



for(i in seq_along(GO_res_list)){
  #for(i in c(7,8)){
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
  
  if(T0_clust < 10){
    T0_clust_idx1 <- GO_matrix2_clust1 == paste0("0",T0_clust)
  } else { T0_clust_idx1 <- GO_matrix2_clust1 == T0_clust}
  T0_clust_idx2 <- GO_matrix2_clust2 == T0_clust
  T0_columns <- c(colnames(GO_matrix2)[T0_clust_idx2],colnames(GO_matrix2)[T0_clust_idx1])
  
  dim(GO_matrix2)
  
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
  rownames(GO_matrix2) <- gsub( "Homo sapiens.*" ,  "" , rownames(GO_matrix2))
  
  
  f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
  ht1 = Heatmap(as.matrix(GO_matrix2),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = T, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = F, cluster_rows = F , row_names_side ="left", width = unit(ncol(GO_matrix2), "cm"),) # width = unit(140, "cm"),
  #ht1
  
  #pdf(paste( Project, "ComplexHeatmap", dataset_name, "Fig__ALL_Pairwise_", GO_res_name, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix2)/2+5), height=nrow(GO_matrix2))
  #par(bg=NA)
  #draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  #dev.off()
  
  dim(GO_matrix2)
  GO_matrix3 <- GO_matrix2[colnames(GO_matrix2) %in% T0_columns]
  dim(GO_matrix3)
  GO_matrix3 <- GO_matrix3[apply(GO_matrix3, 1, function(y) !all(is.na(y))),]
  dim(GO_matrix3)
  GO_matrix3[is.na(GO_matrix3)] <- 1
  
  f1 = colorRamp2( c(0, 0.00001, 0.001, 0.05, 0.5, 1), c("blue4", "darkorchid4", "maroon3", "orchid1", "white", "lightgrey"), space = "RGB") 
  ht1 = Heatmap(as.matrix(GO_matrix3),  col = f1, name = GO_res_name,  row_title = "", column_title = GO_res_name, show_row_names = TRUE, na_col= "lightgrey", heatmap_legend_param = list(title = "Significant terms", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = TRUE, cluster_rows = TRUE , row_names_side ="left", width = unit(ncol(GO_matrix3), "cm")) # width = unit(140, "cm"),
  print(ht1)
  
  pdf(paste( Project, "ComplexHeatmap", "Together", "Fig___T0_to_others_", GO_res_name, "l2fc", l2fc_cutoff, ".pdf", sep="_"), onefile=FALSE, width=(ncol(GO_matrix3)/2+5), height=nrow(GO_matrix3)/3)
  par(bg=NA)
  draw(ht1, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
  
}





message("--------------------------------------------------------------------------------")
message("+                dot plots for selected WP pathways                             ")
message("+-------------------------------------------------------------------------------")

WP_genes_df <- read.csv("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/WP_Gene_list_for_dotplots.csv")

WP_ETC <- na.omit(as.character(unique(WP_genes_df$Electron.Transport.Chain.WP295)))
WP_OxP <- as.character(unique(WP_genes_df$Oxidative.phosphorylation.WP1248))
WP_GaG <- as.character(unique(WP_genes_df$Glycolysis.and.Gluconeogenesis.WP157))
WP_CRP <- as.character(unique(WP_genes_df$Cytoplasmic.Ribosomal.Proteins.WP163))

WP_ETC <- WP_ETC[WP_ETC %in% rownames(GetAssayData(dropseq.integrated, assay = "SCT"))]
WP_OxP <- WP_OxP[WP_OxP %in% rownames(GetAssayData(dropseq.integrated, assay = "SCT"))]
WP_GaG <- WP_GaG[WP_GaG %in% rownames(GetAssayData(dropseq.integrated, assay = "SCT"))]
WP_CRP <- WP_CRP[WP_CRP %in% rownames(GetAssayData(dropseq.integrated, assay = "SCT"))]



pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", Project, "WP_ElectronTransportChain",".pdf", sep="_"), width=length(WP_ETC)/5+2, height=6)
par(bg=NA)
DotPlot(dropseq.integrated, features = WP_ETC , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", Project, "WP_OxidativePhosphorylation",".pdf", sep="_"), width=length(WP_OxP)/5+2, height=6)
par(bg=NA)
DotPlot(dropseq.integrated, features = WP_OxP , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", Project, "WP_GlycolysisAndGlucogenesis",".pdf", sep="_"), width=length(WP_GaG)/5+2, height=6)
par(bg=NA)
DotPlot(dropseq.integrated, features = WP_GaG , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste("DotPlot_matrix.su_realign_Harmony_separate", Project, "WP_CytoplasmicRibosomalProteins",".pdf", sep="_"), width=length(WP_CRP)/5+2, height=6)
par(bg=NA)
DotPlot(dropseq.integrated, features = WP_CRP , cols = c("green", "blue"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()





message("--------------------------------------------------------------------------------")
message("+               Gene cluster semantic similarity measurement                    ")
message("+-------------------------------------------------------------------------------")
# https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html

library(biomaRt)
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene_id'), mart = ensembl)          


library(GOSemSim)
mmGOBP <- godata('org.Mm.eg.db', ont="BP",  keytype = "SYMBOL")
mmGOMF <- godata('org.Mm.eg.db', ont="MF", keytype = "SYMBOL")



Project <- "CTR_dscj1_HARMONY_500g_res0.8_markers"

markerDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g/"
setwd(markerDir)

result_list <- readRDS("Harmony.orig_500g_Markers/result_list_harmony.orig_500g_res0.8.Rds")
names(result_list)
head(result_list[[1]][[3]][[1]],2)



l2fc_cutoff <- 1

extract_genes_from_res_list <- function(x,y){
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
      print(tmp_df$Row.names) 
      return(tmp_df$Row.names)
    }
  }
}


marker_gene_list <- list()
for (z in seq_along(result_list)){
  gene_list <- list()
  res <- result_list[[z]]
  res_name <- names(result_list)[[z]]
  for (w in 1:length(res)){
    gene_list[[w]] <- extract_genes_from_res_list(w, res)
  }
  marker_gene_list[[z]] <- gene_list
}

#marker_gene_list_l2fc0.25 <- marker_gene_list
#marker_gene_list_l2fc0.6 <- marker_gene_list
marker_gene_list_l2fc1 <- marker_gene_list

###marker_genes_l2fc1.5 <- unique(unlist(unlist(marker_gene_list)))


# goSim() function calculates semantic similarity between two GO terms
# mgoSim() function calculates semantic similarity between two sets of GO terms.
# geneSim() to calculate semantic similarity between two gene products, 
# mgeneSim() to calculate semantic similarity among multiple gene products.
# clusterSim() calculating semantic similarity between two gene clusters 
# mclusterSim() function measuring semantic similarity among multiple gene clusters using 


cluster_genes <- function(marker_gene_list_X, first_clust, second_clust) {
  i <- as.numeric(first_clust) + 1
  j <- as.numeric(second_clust) + 2 - i
  if (first_clust > second_clust) { print("First cluster has to be smaller than second!!!") } else {
    print(paste0("first clust coord: ", i, ", second clust coord: ", j))
    clust_genes <- marker_gene_list_X[[i]][[j]]
    return(clust_genes)
  }
}

cluster_genes(marker_gene_list_l2fc0.25, first_clust= 0, second_clust= 21)

gs1 <- cluster_genes(marker_gene_list_l2fc0.25, first_clust= 0, second_clust= 19)
gs2 <- cluster_genes(marker_gene_list_l2fc0.25, first_clust= 0, second_clust= 7)
clusterSim(gs1, gs2, semData=mmGOBP, measure="Wang", combine="BMA") # 0.7  so works!

GO_res_list[[1]][[1]][[1]]

mclusterSim(clusters, semData=hsGO, measure="Wang", combine="BMA")





# russell::: use mgoSim() in an all-against-all pairwise comparison for all UMAP clusters. This was used in the two different heatmaps, but most usefully to get the similarity scores which I put in a matrix and plotted.


GOBP_df <- read.csv("Harmony.orig_500g_Markers_GO/NewAlign51_SCT_HARMONY_500g__res_0.8__GOBP_matrix_l2fc0.25_.csv", row.names = "X")
head(GOBP_df,2)
rownames(GOBP_df) <- GOBP_df[,1]


get_GO_for_clust_pair <- function(GOBP_df_X, first_clust, second_clust){
  if (first_clust > second_clust){ print("First cluster has to be smaller than second!!!") } 
  else {
    #if (first_clust < 10) { 
    #  first_clust <- paste0(0, first_clust)
    #}
    sel_col <- paste0("X", first_clust, ".vs.", second_clust)
    rownames(GOBP_df_X) <- gsub( ".*\\(", "", rownames(GOBP_df_X) )
    rownames(GOBP_df_X) <- gsub( "\\)", "", rownames(GOBP_df_X) )
    idx <- which(colnames(GOBP_df_X) %in% sel_col)
    GO_df_filtered <- na.omit(GOBP_df_X[, c(1, idx)])
    return(rownames(GO_df_filtered))
  }
}

# testing for single pair of GO terms cluster pair::: 

GO_terms_test1 <- get_GO_for_clust_pair(GOBP_df, first_clust = 1, second_clust = 2)
GO_terms_test2 <- get_GO_for_clust_pair(GOBP_df, first_clust = 1, second_clust = 3)
mgoSim(GO_terms_test1, GO_terms_test2, semData=mmGOBP, measure="Wang", combine="BMA")


# set up table to fill in with SCORES:
clusters <- c(0:20) # for together dataset! - cluster 20 may need to be removed !
clusterPairs <- expand.grid(clusters,clusters)
clusterPairs <- clusterPairs[clusterPairs$Var1 != clusterPairs$Var2 , ]
clusterPairs <- clusterPairs[clusterPairs$Var1 < clusterPairs$Var2 , ]
clusterPairs$name <- paste0("X", clusterPairs$Var1, ".vs.", clusterPairs$Var2)
colnames(GOBP_df) <- gsub( "X0",  "X", colnames(GOBP_df))

GO_semSim_tbl <- dplyr::as_tibble(matrix(0, ncol = length(clusterPairs$name), nrow = length(clusterPairs$name) ))
rownames(GO_semSim_tbl) <- clusterPairs$name
colnames(GO_semSim_tbl) <- clusterPairs$name
GO_semSim_tbl <- as.data.frame(GO_semSim_tbl)


calculate_mgoSim_for_2_clustPairs <- function(df, clusterPair1, clusterPair2, semData=mmGOBP){
  # clusterPair1 and clusterPair2 provided as vector or 2 values in increasing order:  e.g. c(1,2)
  GO_terms_clustPair1 <- get_GO_for_clust_pair(GOBP_df, first_clust = clusterPair1[[1]], second_clust = clusterPair1[[2]])
  GO_terms_clustPair2 <- get_GO_for_clust_pair(GOBP_df, first_clust = clusterPair2[[1]], second_clust = clusterPair2[[2]])
  score <- mgoSim(GO_terms_clustPair1, GO_terms_clustPair2, semData=mmGOBP, measure="Wang", combine="BMA")
  pair1Name <- paste0("X", clusterPair1[[1]], ".vs.", clusterPair1[[2]])
  pair2Name <- paste0("X", clusterPair2[[1]], ".vs.", clusterPair2[[2]])
  idx1 <- which(rownames(GO_semSim_tbl) == pair1Name)
  idx2 <- which(colnames(GO_semSim_tbl) == pair2Name)
  GO_semSim_tbl[idx1, idx2] <- score
  print(paste(idx1, "," ,idx2," ", score))
}

calculate_mgoSim_for_2_clustPairs(GOBP_df, clusterPair1= c(1,2), clusterPair2= c(1,3))
GO_semSim_tbl[3, 5]







go1 = c("GO:0004022","GO:0004024","GO:0004174")
go2 = c("GO:0009055","GO:0005515")
mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)
##            GO:0009055 GO:0005515
## GO:0004022      0.368      0.116
## GO:0004024      0.335      0.107
## GO:0004174      0.663      0.119
mgoSim(go1, go2, semData=hsGO, measure="Wang", combine="BMA")
## [1] 0.43







message("--------------------------------------------------------------------------------")
message("+        cluster 3 from T24 compare with together dataset clusters                   ")
message("+-------------------------------------------------------------------------------")

T24_clust3_cells <- readRDS("/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_T24_500g/T24_clust3_cells.Rds")

metaData_together <- dropseq.integrated@meta.data

metaData_together$T24_clust3 <- ifelse(rownames(metaData_together) %in% T24_clust3_cells, "T24_clust3_cells", "other")
colnames(metaData_together)

table(metaData_together[,c(20,24)])




message("+-------------------------------------------------------------------------------")
message("+               Tgfbeta and pathway R I plots                        ")
message("+-------------------------------------------------------------------------------")

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Experiment
markers_T4 <- FindMarkers(dropseq.integrated, ident.1 = "R.4", ident.2= "I.4", logfc.threshold = 0.25, verbose = TRUE)
markers_T1 <- FindMarkers(dropseq.integrated, ident.1 = "R.1", ident.2= "I.1", logfc.threshold = 0.25, verbose = TRUE)
markers_T24 <- FindMarkers(dropseq.integrated, ident.1 = "R.24", ident.2= "I.24", logfc.threshold = 0.25, verbose = TRUE)
markers_T36 <- FindMarkers(dropseq.integrated, ident.1 = "R.36", ident.2= "I.36", logfc.threshold = 0.25, verbose = TRUE)
markers_T48 <- FindMarkers(dropseq.integrated, ident.1 = "R.48", ident.2= "I.48", logfc.threshold = 0.25, verbose = TRUE)

markers_R_I <- FindMarkers(dropseq.integrated, ident.1 = c("R.1", "R.4", "R.36", "R.48"), ident.2= c("I.1", "I.4","I.36", "I.48"), logfc.threshold = 0.25, verbose = TRUE)
markers_R_I <- rownames(markers_R_I)


tbl_markers_freq <- as.data.frame(table(c(rownames(markers_T1), rownames(markers_T4), rownames(markers_T24), rownames(markers_T36), rownames(markers_T48) )) )
tbl_markers_freq <- tbl_markers_freq[order(tbl_markers_freq$Freq, decreasing = TRUE),]
tbl_markers_freq2 <- tbl_markers_freq[tbl_markers_freq$Freq > 2,]

markers_T1$gene <- rownames(markers_T1)
markers_T4$gene <- rownames(markers_T4)
markers_T24$gene <- rownames(markers_T24)
markers_T36$gene <- rownames(markers_T36)
markers_T48$gene <- rownames(markers_T48)
All_R_I_markers <- rbind(markers_T1, markers_T4)
All_R_I_markers <- rbind(All_R_I_markers, markers_T24)
All_R_I_markers <- rbind(All_R_I_markers, markers_T36)
All_R_I_markers <- rbind(All_R_I_markers, markers_T48)
All_R_I_markers <- All_R_I_markers[order(abs(All_R_I_markers$avg_logFC)),]
All_R_I_markers <- All_R_I_markers[!duplicated(All_R_I_markers$gene),]
rownames(All_R_I_markers) <- All_R_I_markers$gene
All_R_I_markers <- All_R_I_markers[All_R_I_markers$gene %in% tbl_markers_freq2$Var1,]
All_R_I_markers <- All_R_I_markers[order((All_R_I_markers$avg_logFC)),]
All_R_I_markers <- All_R_I_markers[All_R_I_markers$gene != "Rpl37",]
All_R_I_markers <- All_R_I_markers[All_R_I_markers$gene != "Rps19",]


pdf(paste("DotPlot_", Project, "DE_at_min3timePoints__l2fc0.25", "_R_vs_I.pdf", sep = "_"), width=9, height=5)
par(bg=NA)
DotPlot(dropseq.integrated, features = All_R_I_markers$gene,   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




#DotPlot(dropseq.integrated, features = tbl_markers_freq[tbl_markers_freq$Freq > 1,]$Var1,   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#DotPlot(dropseq.integrated, features = rownames(markers_R_I),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#DotPlot(dropseq.integrated, features = rownames(R_I_markers),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


markers_T4 <- c("Ncl"  ,   "Ftl1"  ,  "Rpl37" ,  "Klf6" ,   "Serbp1" , "Rpl22" ,  "Krt18" , "Fabp3" ,  "mt-Rnr2")
markers_T1 <- c("Rpl41" , "Rplp2" , "Malat1")
markers_T24 <- c("Ncl" ,   "Rps5" ,  "Id2" ,   "Crip1" , "Malat1" ,"Fabp3") 
markers_R_vs_I <- unique(c(as.character(markers_T1), as.character(markers_T4), as.character(markers_T24)) )
#Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Age
dp1 <- DotPlot(dropseq.integrated, features = markers_R_vs_I,   assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste("DotPlot_", Project, "_", "check_R_vs_I.pdf", sep = "_"), width=9, height=5)
par(bg=NA)
dp1
dev.off()



Idents(dropseq.integrated) <- dropseq.integrated@meta.data[[reso]]

dp2 <- DotPlot(dropseq.integrated, features = c("Tgfb1", "Eomes", "Elf5", "Sox2", "Esrrb", "Cdx2", "Map2k1", "Acvr1" ))

pdf(paste("DotPlot_", Project, "_", "check_R_vs_I_2.pdf", sep = "_"), width=9, height=5)
par(bg=NA)
dp2
dev.off()

dp3 <- DotPlot(dropseq.integrated, features = markers_R_I,   assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste("DotPlot_", Project, "_", "check_R_vs_I_3.pdf", sep = "_"), width=12, height=5)
par(bg=NA)
dp3
dev.off()

DotPlot(dropseq.integrated, features = c("Hand1", "Plet1","Elf5", "Esrrb", "Cdx2", "Eomes", "Tead4", "Sox2", "Yap1", "Krt18",    "Prl8a9", "Prl3d1","Prl2c2", "Pcdh12", "Cdkn1c", "Ets2","Gjb3", "Ascl2", "Tpbpa", "Prdm1",  "Stra13", "Mdfi",      "Tfeb", "Ly6e","Pparg" , "Hif1a",  "Gcm1", "Dlx3", "Ovol2", "Cebpa","Syna", "Junb", "Fosl1", "Arnt") , cols = c("green", "purple"), assay = "SCT", dot.scale = 10, scale = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


Idents(dropseq.integrated) <- paste( "T_", dropseq.integrated@meta.data$Age, sep = "")
Idents(dropseq.integrated) <- factor(Idents(dropseq.integrated) , levels = c( "T_0", "T_1", "T_4", "T_24", "T_36", "T_48") )

plots <- VlnPlot(dropseq.integrated, features = c("Klf6", "Bsg", "Id2", "Crip1", "Serbp1", "Ncl"), split.by = "Treatment", pt.size = 0, combine = FALSE, assay = "SCT")

pdf(paste("DotPlot_", Project, "_", "check_R_vs_I_4.pdf", sep = "_"), width=10, height=16)
par(bg=NA)
wrap_plots(plots = plots, ncol = 1)
dev.off()



#group.by = "celltype", 
DotPlot(dropseq.integrated, features = c("Klf6", "Bsg", "Id2", "Crip1", "Serbp1", "Ncl", "Junb", "Ovol2", "Dlx3",  "Ly6e","Ets2","Gata3", "Ascl2", "Bcor","Tpbpa","Hand1","Hif1a","Gjb3","Tfap2c", "Plac1", "Fgfbp1","Phlda2","Pcdh12","Cdkn1c","Tead4","Tfeb","Stra13","Plet1", "Sparc","Syna", "Gata2","Gcm1", "Tmem37" ,  "Gjb2"  ,   "Car4"  ,   "P4hb"   ,  "Wfdc2" ,   "F2rl1" ,   "Tmem150a", "Igf2" ,    "Peg3" ,    "Igfbp2" ,  "Slc38a4" ),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

table(dropseq.integrated@meta.data[,c("Experiment","Cluster")])


markers_c1 <- FindMarkers(dropseq.integrated, ident.1 = "1", logfc.threshold = 0.4, verbose = TRUE, only.pos = TRUE)
markers_c6 <- FindMarkers(dropseq.integrated, ident.1 = "6", logfc.threshold = 0.4, verbose = TRUE, only.pos = TRUE)
markers_c4 <- FindMarkers(dropseq.integrated, ident.1 = "4", logfc.threshold = 0.4, verbose = TRUE, only.pos = TRUE)
markers_c7 <- FindMarkers(dropseq.integrated, ident.1 = "7", logfc.threshold = 0.25, verbose = TRUE, only.pos = TRUE)
markers_c16 <- FindMarkers(dropseq.integrated, ident.1 = "16", logfc.threshold = 0.4, verbose = TRUE, only.pos = TRUE)
markers_c17 <- FindMarkers(dropseq.integrated, ident.1 = "17", logfc.threshold = 0.4, verbose = TRUE, only.pos = TRUE)
markers_c3_vs_c12 <- FindMarkers(dropseq.integrated, ident.1 = "3", ident.2="12", logfc.threshold = 0.4, only.pos = TRUE)
markers_c2_vs_c12 <- FindMarkers(dropseq.integrated, ident.1 = "2", ident.2="12", logfc.threshold = 0.4, only.pos = TRUE)

# c6 
DotPlot(dropseq.integrated, features = c( "Gata2","Gcm1", "Tmem37" ,  "Gjb2"  ,   "Car4"  ,   "P4hb"   ,  "Wfdc2" ,   "F2rl1" ,   "Tmem150a", "Igf2" ,    "Peg3" ,    "Igfbp2" ,  "Slc38a4", "Fabp3"  ,  "Bsg" , "Maged1"  , "Lgals1"   ,  "Peg10",   "Xist"  ,   "Ldha"   ),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# c1
DotPlot(dropseq.integrated, features = c( "Ascl2","Gjb3","Tfap2c", "Plac1", "Fgfbp1","Phlda2","Hbegf","Tbrg1", "Rbbp7","Hif1a", "Hat1"),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# c4 
DotPlot(dropseq.integrated, features = c("Plac1", "Tfeb", "Ets2", "Ovol2", "Dlx3","Sct"  ,  "Cdkn1c", "Krt19"  ,"Lgals3" ,"Sin3b" , "Krt18" , "Krt8"  , "Tfrc"  , "Ctsb" ,  "Hbegf" , "Psap"  , "Sept4",  "Anxa1", "Las1l" ,   "H2-D1",    "Ralb" ,"Serpine2", "Acadl" ,   "Rsad2" ,   "Hspb1" ,   "Nrk" ,     "Atxn10" ,  "Peg3" ,    "Lgals9",   "Slc16a1",  "Lgals1",   "Bsg" ,     "Cnn2"),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# c16
DotPlot(dropseq.integrated, features = c("Saa3",   "Fabp3" , "Gjb2",   "F2rl1",  "Lgals1", "Igf2" ,  "Psap" ,  "Crip1" , "P4hb"  , "Peg3"  , "Tmem37" ,"Car4"  , "Cdkn1c"),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# c17
DotPlot(dropseq.integrated, features = c("Ccl2",   "Pfn1" , "Mrpl33",   "Fknp3",  "Dynll1", "mt-Cytb" ),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# early diff - clust 2,3 vs 12 T0

DotPlot(dropseq.integrated, features =c("mt-Cytb", "Calm1", "Bsg", "Lgals1", "Phlda2","Rhox6" , "Igf2" ,  "Cdkn1c" ,"Krt18" ,"Tbrg1",  "Malat1", "H2-D1",  "H2-K1" , "Trip12" ,"Xist",   "Anxa2",  "Klhl13", "Rbbp7",  "Krt8" ,  "Ctnna1", "Hat1" ,  "Fabp3"),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))











message("+-------------------------------------------------------------------------------")
message("+       expression of I vs R    plots       ")
message("+-------------------------------------------------------------------------------")


#  all cell expr line plot R and I 
Metadata <- dropseq.integrated@meta.data
Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Experiment


ExprSCT <- GetAssayData(dropseq.integrated, assay = "SCT")
ExprSCT[1:5,1:5]


avg_expt2 = AverageExpression(dropseq.integrated, return.seurat = FALSE) 
avg_expr_expt <- avg_expt2[[2]]

#gene <- "Klf6"


col_expt <- c("0.0"="black", "I.1"="yellow", "R.1"="yellow2", "I.4"="orangered", "R.4"="red", "I.24"="lightblue2", "R.24"="lightblue4", "I.36"="blue", "R.36"="darkblue","I.48"="purple", "R.48"="orchid")  







#  bar plot showing expression through time R next to I

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Age
plots <- VlnPlot(dropseq.integrated, features = c("Id2", "Plac1", "Klf6"), split.by = "Treatment", group.by = "Age",  pt.size = 0, combine = FALSE, slot = "data")
wrap_plots(plots = plots, ncol = 1)

VlnPlot(dropseq.integrated, features = c( "Klf6"), split.by = "Treatment", group.by = "Age", 
        pt.size = 0, combine = FALSE, slot = "data")









message("+     function summarySE      ")


# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     quant = quantile(xx[[col]], na.rm=na.rm, probs = c(0.1, 0.25, 0.75, 0.9))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac)[colnames(datac) == "mean"] <-  measurevar
  #datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  datac$CI_upper <- datac$gene + 0.95*datac$sd/sqrt(datac$N )
  datac$CI_lower <- datac$gene - 0.95*datac$sd/sqrt(datac$N )
  
  
  return(datac)
}
tgc <- summarySE(df_gene, measurevar="gene", groupvars=c("time", "treatment"))





message("+      function:  plt_lin_with_errorbars       ")

plt_lin_with_errorbars <- function(gene, error_type = sd){
  
  df_gene <- as.data.frame(ExprSCT[gene,])
  colnames(df_gene)[1] <- "gene"
  
  df_gene$cell <- rownames(df_gene)
  df_gene$Experiment <- Metadata[match(rownames(df_gene), rownames(Metadata)),]$Experiment
  df_gene$Experiment <- as.character(df_gene$Experiment )
  df_gene$Experiment[df_gene$Experiment == "0.0"] <- "R.0"
  df_gene_dup_T0 <- df_gene[df_gene$Experiment == "R.0",]
  df_gene_dup_T0$Experiment[df_gene_dup_T0$Experiment == "R.0"] <- "I.0"
  
  df_gene <- rbind(df_gene, df_gene_dup_T0)
  
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  df_gene$time <- as.numeric(as.character(df_gene$time))
  
  head(df_gene)
  df_gene <- df_gene[order(df_gene$time),]
  
  tgc <- summarySE(df_gene, measurevar="gene", groupvars=c("time", "treatment"))
  
  
  # Standard error of the mean
  plt <- ggplot(tgc, aes(x=jitter(time), y=gene, colour=treatment)) + 
    #geom_errorbar(aes(ymin=gene-ci, ymax=gene+ci), width=.1) +
    geom_line() + theme_classic() +  scale_color_manual(values=c("R"= 'red', "I"= 'blue')) + 
    geom_point() + geom_errorbar(aes(ymin=gene-tgc[[error_type]], ymax=gene+tgc[[error_type]]), width=.2,  position=position_dodge(.9)) +
    ylab("Average Gene Expression") + xlab("Time (h)") + ggtitle(gene) + geom_ribbon(aes(ymin=gene-ci, ymax=gene+ci), linetype=2, alpha=0.1)
  
  return(plt)
}

plt_lin_with_errorbars("Klf6", error_type = "sd")
plt_lin_with_errorbars("Klf6", error_type = "se")
plt_lin_with_errorbars("Klf6", error_type = "ci")
plt_lin_with_errorbars("Plac1")
plt_lin_with_SD("Id2")




message("+      function:  plt_lin_with_quantiles       ")

plt_lin_with_quantiles <- function(gene){
  
  df_gene <- as.data.frame(ExprSCT[gene,])
  colnames(df_gene)[1] <- "gene"
  
  df_gene$cell <- rownames(df_gene)
  df_gene$Experiment <- Metadata[match(rownames(df_gene), rownames(Metadata)),]$Experiment
  df_gene$Experiment <- as.character(df_gene$Experiment )
  df_gene$Experiment[df_gene$Experiment == "0.0"] <- "R.0"
  df_gene_dup_T0 <- df_gene[df_gene$Experiment == "R.0",]
  df_gene_dup_T0$Experiment[df_gene_dup_T0$Experiment == "R.0"] <- "I.0"
  
  df_gene <- rbind(df_gene, df_gene_dup_T0)
  
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  df_gene$time <- as.numeric(as.character(df_gene$time))
  
  head(df_gene)
  df_gene <- df_gene[order(df_gene$time),]
  
  tgc <- summarySE(df_gene, measurevar="gene", groupvars=c("time", "treatment"))
  
  
  # Standard error of the mean
  plt <- ggplot(tgc, aes(x=jitter(time), y=gene, colour=treatment)) + 
    #geom_errorbar(aes(ymin=gene-ci, ymax=gene+ci), width=.1) +
    geom_line() + theme_classic() +  scale_color_manual(values=c("R"= 'red', "I"= 'blue')) + 
    geom_point() + geom_errorbar(aes(ymin=tgc[[7]], ymax=tgc[[8]]), width=.2,  position=position_dodge(.9)) +
    ylab("Average Gene Expression") + xlab("Time (h)") + ggtitle(gene)
  #+ geom_ribbon(aes(ymin=gene-sd, ymax=gene+sd), linetype=2, alpha=0.1)
  
  return(plt)
}

plt_lin_with_quantiles("Klf6")
plt_lin_with_quantiles("Id2")
plt_lin_with_quantiles("Plac1")




message("+      function:  plt_ridges       ")

plt_ridges <- function(gene){
  
  df_gene <- as.data.frame(ExprSCT[gene,])
  colnames(df_gene)[1] <- "gene"
  
  df_gene$cell <- rownames(df_gene)
  df_gene$Experiment <- Metadata[match(rownames(df_gene), rownames(Metadata)),]$Experiment
  df_gene$Experiment <- as.character(df_gene$Experiment )
  df_gene$Experiment[df_gene$Experiment == "0.0"] <- "R.0"
  df_gene_dup_T0 <- df_gene[df_gene$Experiment == "R.0",]
  df_gene_dup_T0$Experiment[df_gene_dup_T0$Experiment == "R.0"] <- "I.0"
  
  df_gene <- rbind(df_gene, df_gene_dup_T0)
  
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  df_gene$time <- as.numeric(as.character(df_gene$time))
  
  head(df_gene)
  df_gene <- df_gene[order(df_gene$time),]
  df_gene$time <- factor(df_gene$time, levels = c("48", "36", "24", "4", "1", "0"))
  plt_ridges <- ggplot(df_gene, aes(x = gene, y = time, fill = treatment)) + geom_density_ridges2( alpha = 0.4) + scale_fill_manual(values =c("R" = "red", "I" = "blue")) +  xlab("Gene Expression") + ylab("Time (h)") + ggtitle(gene)
  return(plt_ridges)
}

plt_ridges("Klf6")
plt_ridges("Plac1")
plt_ridges("Id2")



plt_ridges2 <- function(gene){
  
  df_gene <- as.data.frame(ExprSCT[gene,])
  colnames(df_gene)[1] <- "gene"
  
  df_gene$cell <- rownames(df_gene)
  df_gene$Experiment <- Metadata[match(rownames(df_gene), rownames(Metadata)),]$Experiment
  df_gene$Experiment <- as.character(df_gene$Experiment )
  df_gene$Experiment[df_gene$Experiment == "0.0"] <- "R.0"
  df_gene_dup_T0 <- df_gene[df_gene$Experiment == "R.0",]
  df_gene_dup_T0$Experiment[df_gene_dup_T0$Experiment == "R.0"] <- "I.0"
  
  df_gene <- rbind(df_gene, df_gene_dup_T0)
  
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  df_gene$time <- as.numeric(as.character(df_gene$time))
  
  head(df_gene)
  df_gene <- df_gene[order(df_gene$time),]
  df_gene$time <- factor(df_gene$time, levels = c("0", "1", "4", "24", "36", "48"))
  plt_ridges <- ggplot(df_gene, aes(x = gene, y = time, fill = treatment)) + geom_density_ridges2( alpha = 0.4) + scale_fill_manual(values =c("R" = "red", "I" = "blue")) +  xlab("Gene Expression") + ylab("Time (h)") + ggtitle(gene) + coord_flip() + theme(legend.title = element_blank()) + theme(legend.position = "none")
  
  return(plt_ridges)
}

pdf(paste("PlotRidges", Project, "LEGEND",".pdf", sep="_"), width=6, height=4)
par(bg=NA)
plt_ridges2("Klf6")
dev.off()





ggplot(data=df_gene, aes(x=time, y=jitter(gene), colour = treatment)) +
  geom_line(linetype = "solid") + geom_point() + theme_classic() +   
  scale_color_manual(values=c("R"= 'red', "I"= 'blue')) + 
  ylab("Average Gene Expression") + xlab("Time (h)") + ggtitle(gene) #+ geom_quantile(quantiles = c(0.25,0.75)) #+ geom_smooth(se=TRUE)
# + geom_ribbon(aes(ymin=df_gene$lower, ymax=df_gene$upper), linetype=2, alpha=0.1)






message("+      function:  gene_linear_plot       ")

#  avg expr line plot R and I 

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Experiment
avg_expt = AverageExpression(dropseq.integrated, return.seurat = TRUE) # , add.ident = 'Experiment'
avg_expt2 = AverageExpression(dropseq.integrated, return.seurat = FALSE, slot = "data", assays = "SCT" ) 

avg_expr_expt <- avg_expt2[[2]]


gene_linear_plot <- function(gene){
  
  df_gene <- avg_expr_expt[gene,]
  colnames(df_gene)[colnames(df_gene) == "0.0"] <- "R.0"
  df_gene$`I.0` <- df_gene$`R.0`
  
  df_gene <- as.data.frame(t(df_gene))
  df_gene$time <- gsub( ".*\\.", "", rownames(df_gene))
  df_gene$treatment <- gsub( "\\..*", "", rownames(df_gene))
  df_gene$time <- as.numeric(as.character(df_gene$time))
  colnames(df_gene)[1] <- "gene"
  
  plt_lin <- ggplot(data=df_gene, aes(x=time, y=gene, colour = treatment)) +
    geom_line(linetype = "solid") + geom_point() + theme_classic() +   
    scale_color_manual(values=c("R"= 'red', "I"= 'blue')) + 
    ylab("Avg. Gene Expression") + xlab("Time (h)") + ggtitle(gene) + 
    theme(legend.title = element_blank()) + theme(legend.position = "none")
  return(plt_lin)
}

gene_linear_plot("Id2")
gene_linear_plot("Klf6")
gene_linear_plot("Crip1")
gene_linear_plot("Bsg")
gene_linear_plot("Ncl")
gene_linear_plot("Serbp1")

pdf(paste("gene_linear_plot", Project, "LEGEND",".pdf", sep="_"), width=6, height=3)
par(bg=NA)
gene_linear_plot("Id2")
dev.off()



gene_linear_plot2 <- function(gene){
  
  df_gene <- as.data.frame(ExprSCT[gene,])
  colnames(df_gene)[1] <- "gene"
  
  df_gene$cell <- rownames(df_gene)
  df_gene$Experiment <- Metadata[match(rownames(df_gene), rownames(Metadata)),]$Experiment
  df_gene$Experiment <- as.character(df_gene$Experiment )
  df_gene$Experiment[df_gene$Experiment == "0.0"] <- "R.0"
  df_gene_dup_T0 <- df_gene[df_gene$Experiment == "R.0",]
  df_gene_dup_T0$Experiment[df_gene_dup_T0$Experiment == "R.0"] <- "I.0"
  
  df_gene <- rbind(df_gene, df_gene_dup_T0)
  
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  df_gene$time <- as.numeric(as.character(df_gene$time))
  
  head(df_gene)
  df_gene <- df_gene[order(df_gene$time),]
  df_gene$time <- factor(df_gene$time, levels = c("0", "1", "4", "24", "36", "48"))
  
  mean_expr <- aggregate(df_gene[, 1], list(df_gene$Experiment), mean)
  rownames(mean_expr) <- mean_expr[,1]
  mean_expr$treatment <- gsub( "\\..*", "", mean_expr$Group.1)
  mean_expr$time <- gsub( ".*\\.", "", mean_expr$Group.1)
  mean_expr$time <- as.numeric(as.character((mean_expr$time)))
  
  plt_lin <- ggplot(data=mean_expr, aes(x=time, y=x, colour = treatment)) +
    geom_line(linetype = "solid") + geom_point() + theme_classic() +   
    scale_color_manual(values=c("R"= 'red', "I"= 'blue')) + 
    ylab("Avg. Gene Expression") + xlab("Time (h)") + ggtitle(gene) + 
    theme(legend.title = element_blank()) + theme(legend.position = "none")
  return(plt_lin)
}
gene_linear_plot2("Klf6")








ExprSCT <- GetAssayData(dropseq.integrated, assay = "SCT", slot = "data")
#Metadata <- dropseq.integrated@meta.data

avg_expt2 = AverageExpression(dropseq.integrated, return.seurat = FALSE, slot = "data", assays = "SCT" ) 
avg_expr_expt <- avg_expt2[[1]]


plt_ridges_with_linear <- function(gene){
  
  #ExprSCT <- GetAssayData(dropseq.integrated, assay = "SCT", slot = "data")
  #Metadata <- dropseq.integrated@meta.data

  #avg_expt2 = AverageExpression(dropseq.integrated, return.seurat = FALSE, slot = "data", assays = "SCT" ) 
  #avg_expr_expt <- avg_expt2[[1]]
  
  df_gene <- as.data.frame(ExprSCT[gene,])
  colnames(df_gene)[1] <- "gene"
  
  df_gene$cell <- rownames(df_gene)
  df_gene$Experiment <- Metadata[match(rownames(df_gene), rownames(Metadata)),]$Experiment
  df_gene$Experiment <- as.character(df_gene$Experiment )
  df_gene$Experiment[df_gene$Experiment == "0.0"] <- "R.0"
  df_gene_dup_T0 <- df_gene[df_gene$Experiment == "R.0",]
  df_gene_dup_T0$Experiment[df_gene_dup_T0$Experiment == "R.0"] <- "I.0"
  
  df_gene <- rbind(df_gene, df_gene_dup_T0)
  
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  df_gene$time <- as.numeric(as.character(df_gene$time))
  
  head(df_gene)
  df_gene <- df_gene[order(df_gene$time),]
  df_gene$time <- factor(df_gene$time, levels = c("0", "1", "4", "24", "36", "48"))
  
  mean_gene <- avg_expr_expt[gene,]
  colnames(mean_gene)[colnames(mean_gene) == "0.0"] <- "R.0"
  mean_gene$`I.0` <- mean_gene$`R.0`
  
  mean_gene <- as.data.frame(t(mean_gene))
  mean_gene$time <- gsub( ".*\\.", "", rownames(mean_gene))
  mean_gene$treatment <- gsub( "\\..*", "", rownames(mean_gene))
  mean_gene$time <- as.numeric(as.character(mean_gene$time))
  colnames(mean_gene)[1] <- "gene"
  #mean_gene$time <- as.numeric(as.character(mean_gene$time))
  mean_gene$time <- as.character(mean_gene$time)
  mean_gene$time <- factor(mean_gene$time, levels = c("0", "1", "4", "24", "36", "48"))
  mean_gene <- mean_gene[order(mean_gene$time),]
  
  df_gene$mean_expression <- mean_gene[match(df_gene$Experiment, rownames(mean_gene)),]$gene

  mean_expr <- aggregate(df_gene[, 1], list(df_gene$Experiment), mean)
  mean_gene$mean_expression <- mean_expr[match(rownames(mean_gene), mean_expr$Group.1),]$x
  
  plt_ridges <- ggplot(df_gene, aes(x = gene, y = time, fill = treatment)) + geom_density_ridges2( alpha = 0.4, quantile_lines=TRUE, quantile_fun=function(x,...)mean(x)) + scale_fill_manual(values =c("R" = "red", "I" = "blue")) +  xlab("Gene Expression") + ylab("Time (h)") + ggtitle(gene) + coord_flip() + theme(legend.title = element_blank()) + theme(legend.position = "none") + 
    geom_path(data=mean_gene, mapping=aes(x = mean_expression, y = time, group = treatment, colour = treatment)) +
    scale_color_manual(values =c("R" = "red", "I" = "blue")) 
  return(plt_ridges)
}


#pdf(paste("PlotRidges", Project, "LEGEND",".pdf", sep="_"), width=6, height=4)
#par(bg=NA)
plt_ridges_with_linear("Klf6")
#dev.off()






FeaturePlot(dropseq.integrated, features = c("Gcm1", "E2f8","Atrx"))

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Cluster
DotPlot(dropseq.integrated, features =c("Ascl2","Tfap2c","Arid3a",  "Hic2","Hif1a",     "Bsg","E2f8","Gcm1","Msx2", "Crip1",     "Cdx2","Tead4","Elf5","Elf2","Atrx","Id2","Ncl","Serbp1","Mafk","Tsc22d1","Krt18", "Ets2",  "Bhlhe40","Fbxo21","Creb3l2","Abca4","Ctcf", "cFos", "Irf2", "Maff", "Mef2d", "Meis1", "Pou3f1", "Cbfa2t3", "Foxj2", "Tbx20",  "Hopx", "Smad6", "Klf7","Klf6","Lrrfip1","Bbx", "Pcgf5"),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



pdf(paste("gene_linear_plot", Project, "R_vs_I_genes",".pdf", sep="_"), width=10, height=5)
par(bg=NA)
plot_grid(gene_linear_plot("Id2"), gene_linear_plot("Ncl"), gene_linear_plot("Serbp1"), 
          gene_linear_plot("Klf6"),gene_linear_plot("Crip1"), gene_linear_plot("Bsg"),  ncol = 3)
dev.off()





#  avg expr box plot R and I (per sample for stdev)

sampleTable_red <- unique(dropseq.integrated@meta.data[,c("orig.ident", "Experiment")])

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$orig.ident
avg_sample = AverageExpression(dropseq.integrated, return.seurat = TRUE) 
avg_sample2 = AverageExpression(dropseq.integrated, return.seurat = FALSE)

avg_expr_sample <- avg_sample2[[2]]




message("+      function:  gene_box_plot       ")

gene_box_plot_per_sample <- function(gene) {
  
  df_gene <- avg_expr_sample[gene,]
  
  df_gene <- as.data.frame(t(df_gene))
  df_gene$Experiment <- sampleTable_red[match(rownames(df_gene), sampleTable_red$orig.ident),]$Experiment
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  colnames(df_gene)[1] <- "gene"
  df_gene$time <- factor(df_gene$time, levels = c(0,1,4,24,36,48))
  
  plt_box <- ggplot(data=df_gene, aes(x=time, y=gene, fill = treatment)) +
    geom_boxplot()  + theme_classic() +   geom_line() +
    scale_fill_manual(values=c("R"= 'red', "I"= 'blue', "0" = "green")) + 
    ylab("Average Gene Expression") + xlab("Time (h)") + ggtitle(gene)
  return(plt_box)
}


gene_box_plot <- function(gene) {
  
  df_gene <- as.data.frame(ExprSCT[gene,])
  colnames(df_gene)[1] <- "gene"
  
  df_gene$cell <- rownames(df_gene)
  df_gene$Experiment <- Metadata[match(rownames(df_gene), rownames(Metadata)),]$Experiment
  df_gene$Experiment <- as.character(df_gene$Experiment )
  df_gene$Experiment[df_gene$Experiment == "0.0"] <- "R.0"
  df_gene_dup_T0 <- df_gene[df_gene$Experiment == "R.0",]
  df_gene_dup_T0$Experiment[df_gene_dup_T0$Experiment == "R.0"] <- "I.0"
  
  df_gene <- rbind(df_gene, df_gene_dup_T0)
  
  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
  df_gene$time <- as.numeric(as.character(df_gene$time))
  
  head(df_gene)
  df_gene <- df_gene[order(df_gene$time),]
#  df_gene$time <- factor(df_gene$time, levels = c("48", "36", "24", "4", "1", "0"))
  
#  df_gene <- as.data.frame(t(df_gene))
#  df_gene$Experiment <- sampleTable_red[match(rownames(df_gene), sampleTable_red$orig.ident),]$Experiment
#  df_gene$time <- gsub( ".*\\.", "", df_gene$Experiment)
#  df_gene$treatment <- gsub( "\\..*", "", df_gene$Experiment)
#  colnames(df_gene)[1] <- "gene"
  df_gene$time <- factor(df_gene$time, levels = c(0,1,4,24,36,48))
  
  plt_box <- ggplot(data=df_gene, aes(x=time, y=gene, fill = treatment)) +
    geom_boxplot()  + theme_classic() +   geom_line() +
    scale_fill_manual(values=c("R"= 'red', "I"= 'blue', "0" = "green")) + 
    ylab("Average Gene Expression") + xlab("Time (h)") + ggtitle(gene)
  return(plt_box)
}

gene_box_plot("Id2")
gene_box_plot("Klf6")
gene_box_plot("Crip1")
gene_box_plot("Bsg")
gene_box_plot("Ncl")
gene_box_plot("Serbp1")


pdf(paste("gene_box_plot", Project, "R_vs_I_genes",".pdf", sep="_"), width=10, height=5)
par(bg=NA)
plot_grid(gene_box_plot("Id2"), gene_box_plot("Ncl"), gene_box_plot("Serbp1"), 
          gene_box_plot("Klf6"), gene_box_plot("Crip1"), gene_box_plot("Bsg"),  ncol = 3)
dev.off()



# avg expr I / R



gene_I_divby_R_plot <- function(gene) {
  
  df_gene <- avg_expr_expt[gene,]
  colnames(df_gene)[colnames(df_gene) == "0.0"] <- "R.0"
  df_gene$`I.0` <- df_gene$`R.0`
  
  df_gene <- as.data.frame(t(df_gene))
  df_gene$time <- gsub( ".*\\.", "", rownames(df_gene))
  df_gene$treatment <- gsub( "\\..*", "", rownames(df_gene))
  df_gene$time <- as.numeric(as.character(df_gene$time))
  colnames(df_gene)[1] <- "gene"
  #df_gene$Experiment <- rownames(df_gene)
  #df_gene_melt <- reshape2::melt(df_gene, id.vars = c( "time","gene"))
  df_I <- df_gene[df_gene$treatment == "I",]
  df_R <- df_gene[df_gene$treatment == "R",]
  df_I <- df_I[order(df_I$time),]
  df_R <- df_R[order(df_R$time),]
  colnames(df_I)[1] <- "gene_I"
  colnames(df_R)[1] <- "gene_R"
  df_plt <- cbind(df_I[,c(1,2)], df_R[,c(1,3)])
  df_plt$gene_ratio <- df_plt$gene_I / df_plt$gene_R
  
  plt_ratio_lin <- ggplot(data=df_plt, aes(x=time, y=gene_ratio, colour = treatment)) +
    geom_line(linetype = "solid") + geom_point() + theme_classic() +   
    scale_color_manual(values=c("green")) + theme(legend.position = "none") +
    ylab("Gene Expression Ratio (I / R)") + xlab("Time (h)") + ggtitle(gene) + 
    geom_hline(yintercept= 1) + ylim( 1- (max(df_plt$gene_ratio) ), max(df_plt$gene_ratio) + 0.3)
  return(plt_ratio_lin)
}

gene_I_divby_R_plot("Id2")
gene_I_divby_R_plot("Klf6")
gene_I_divby_R_plot("Crip1")
gene_I_divby_R_plot("Bsg")
gene_I_divby_R_plot("Ncl")
gene_I_divby_R_plot("Serbp1")


pdf(paste("gene_I_divby_R_plot", Project, "R_vs_I_genes",".pdf", sep="_"), width=10, height=5)
par(bg=NA)
plot_grid(gene_I_divby_R_plot("Id2"), gene_I_divby_R_plot("Ncl"), gene_I_divby_R_plot("Serbp1"), 
          gene_I_divby_R_plot("Klf6"),gene_I_divby_R_plot("Crip1"), gene_I_divby_R_plot("Bsg"),  ncol = 3)
dev.off()




Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Treatment
R_I_markers <- FindMarkers(dropseq.integrated, ident.1 = "I", ident.2= "R", logfc.threshold = 0.25, verbose = TRUE)

Idents(dropseq.integrated) <- dropseq.integrated@meta.data$Experiment
DotPlot(dropseq.integrated, features =rownames(R_I_markers),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste("DotPlot", Project, "R_I_markers", "l2fc0.25", ".pdf", sep="_"), width=12, height=4.5)
par(bg=NA)
DotPlot(dropseq.integrated, features =rownames(R_I_markers),   assay = "SCT", dot.scale = 10, scale = TRUE, cols = c("green", "purple")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


pdf(paste("gene_linear_plot", Project, "R_I_markers",".pdf", sep="_"), width=25, height=8)
par(bg=NA)
plot_grid( gene_linear_plot("Id2"), gene_linear_plot("Ncl"), gene_linear_plot("Serbp1"),  gene_linear_plot("Phlda2"), 
           gene_linear_plot("Top2a"), gene_linear_plot("Npm1"), gene_linear_plot("Gas5"),
           gene_linear_plot("Klf6"),gene_linear_plot("Crip1"), gene_linear_plot("Bsg"), gene_linear_plot("Dppa1"), 
           gene_linear_plot("Serpine2"), gene_linear_plot("B2m"),  gene_linear_plot("Rsad2"),  
           gene_linear_plot("Lgals3"), gene_linear_plot("Acadl"), gene_linear_plot("H2-K1"), gene_linear_plot("Psap"),
           gene_linear_plot("Krt8"), gene_linear_plot("Ctr9"), gene_linear_plot("Ralb"), ncol = 7)
dev.off()


#"Fabp3"    "Psap"     "Lgals1"   "Dppa1"    "Lgals3"   "Krt18"    "Krt8"     "Sct"      "Cdkn1c"   "Ralb"     "Peg3"     "H2-D1"   
# "Klf6"     "Ctsb"     "Acadl"    "Crip1"    "P4hb"     "Anxa2"    "Serpine2" "H2-K1"    "Igf2"     "Bsg"      "Krt19"    "Rsad2"   
# "Gjb2"     "Ctr9"     "Malat1"   "B2m"      "L1td1"    "Dstn"     "Hspb1"    "Ctnna1"   "Slc2a1"   "Cnn2"     "Gas5"     "Top2a"   
# "Ranbp1"   "Npm1"     "Id2"      "Ncl"      "Phlda2" 










# testing performance of different plot types:::

plot_grid(plt_lin_with_errorbars("Id2", error_type = "sd"), plt_lin_with_errorbars("Ncl", error_type = "sd"), plt_lin_with_errorbars("Serbp1", error_type = "sd"), 
          plt_lin_with_errorbars("Klf6", error_type = "sd"), plt_lin_with_errorbars("Crip1", error_type = "sd"), plt_lin_with_errorbars("Bsg", error_type = "sd"),  ncol = 3)


plot_grid(plt_lin_with_errorbars("Id2", error_type = "ci"), plt_lin_with_errorbars("Ncl", error_type = "ci"), plt_lin_with_errorbars("Serbp1", error_type = "ci"), 
          plt_lin_with_errorbars("Klf6", error_type = "ci"), plt_lin_with_errorbars("Crip1", error_type = "ci"), plt_lin_with_errorbars("Bsg", error_type = "ci"),  ncol = 3)


plot_grid(gene_linear_plot("Id2"), gene_linear_plot("Ncl"), gene_linear_plot("Serbp1"), 
          gene_linear_plot("Klf6"), gene_linear_plot("Crip1"), gene_linear_plot("Bsg"),  ncol = 3)

plot_grid(plt_ridges("Id2"), plt_ridges("Ncl"), plt_ridges("Serbp1"), 
          plt_ridges("Klf6"), plt_ridges("Crip1"), plt_ridges("Bsg"),  ncol = 3)

plot_grid(plt_ridges2("Id2"), plt_ridges2("Ncl"), plt_ridges2("Serbp1"), 
          plt_ridges2("Klf6"), plt_ridges2("Crip1"), plt_ridges2("Bsg"),  ncol = 3)





plot_grid(gene_box_plot("Id2"), gene_box_plot("Ncl"), gene_box_plot("Serbp1"), 
          gene_box_plot("Klf6"), gene_box_plot("Crip1"), gene_box_plot("Bsg"),  ncol = 3)

plot_grid(gene_box_plot_per_sample("Id2"), gene_box_plot_per_sample("Ncl"), gene_box_plot_per_sample("Serbp1"), 
          gene_box_plot_per_sample("Klf6"), gene_box_plot_per_sample("Crip1"), gene_box_plot_per_sample("Bsg"),  ncol = 3)

plot_grid(gene_I_divby_R_plot("Id2"), gene_I_divby_R_plot("Ncl"), gene_I_divby_R_plot("Serbp1"), 
          gene_I_divby_R_plot("Klf6"),gene_I_divby_R_plot("Crip1"), gene_I_divby_R_plot("Bsg"),  ncol = 3)

plot_grid(plt_lin_with_quantiles("Id2"), plt_lin_with_quantiles("Ncl"), plt_lin_with_quantiles("Serbp1"), 
          plt_lin_with_quantiles("Klf6"),plt_lin_with_quantiles("Crip1"), plt_lin_with_quantiles("Bsg"),  ncol = 3)





plt_lin_with_errorbars("Klf6", error_type = "sd")
plt_lin_with_errorbars("Klf6", error_type = "se")
plt_lin_with_errorbars("Klf6", error_type = "ci")
plt_lin_with_errorbars("Plac1")
plt_lin_with_SD("Id2")



gene_linear_plot("Id2")
gene_linear_plot("Klf6")
gene_linear_plot("Crip1")
gene_linear_plot("Bsg")
gene_linear_plot("Ncl")
gene_linear_plot("Serbp1")

gene_linear_plot("Cdx2")
gene_linear_plot("Tead4")
gene_linear_plot("Elf5")
gene_linear_plot("Krt18")
gene_linear_plot("Tfap2c")
gene_linear_plot("Msx2")
gene_linear_plot("Atrx")
gene_linear_plot("Mafk")
gene_linear_plot("Elf2")
gene_linear_plot("Bhlhe40")
gene_linear_plot("Tsc22d1")
gene_linear_plot("Klf7")
gene_linear_plot("Creb3l2")
gene_linear_plot("")
gene_linear_plot("Abca4")
gene_linear_plot("Fbxo21")
gene_linear_plot("Gcm1")
gene_linear_plot("E2f8")
gene_linear_plot("Atrx")


plt_ridges("Klf6")
plt_ridges("Plac1")
plt_ridges("Id2")



All_R_I_markers$gene

pdf(paste("GeneLinearPlot", Project, "16x_top_R_I_markers", "l2fc0.25", ".pdf", sep="_"), width=10, height=10)
par(bg=NA)
plot_grid(gene_linear_plot("Phlda2"), gene_linear_plot("Ncl"), gene_linear_plot("Npm1"), gene_linear_plot("Id2"), gene_linear_plot("Paics"), 
          gene_linear_plot("G3bp2"), gene_linear_plot("Rpl4"), gene_linear_plot("Cenpf"), gene_linear_plot("Klf6"), gene_linear_plot("B2m"),
          gene_linear_plot("Serpine2"), gene_linear_plot("Krt18"), gene_linear_plot("Acadl"), gene_linear_plot("Lgals1"), gene_linear_plot("Ralb"), 
          gene_linear_plot("Dppa1"),  ncol = 4)
dev.off()

pdf(paste("GeneLinearPlot2", Project, "16x_top_R_I_markers", "l2fc0.25", ".pdf", sep="_"), width=10, height=10)
par(bg=NA)
plot_grid(gene_linear_plot2("Phlda2"), gene_linear_plot2("Ncl"), gene_linear_plot2("Npm1"), gene_linear_plot2("Id2"), gene_linear_plot2("Paics"), 
          gene_linear_plot2("G3bp2"), gene_linear_plot2("Rpl4"), gene_linear_plot2("Cenpf"), gene_linear_plot2("Klf6"), gene_linear_plot2("B2m"),
          gene_linear_plot2("Serpine2"), gene_linear_plot2("Krt18"), gene_linear_plot2("Acadl"), gene_linear_plot2("Lgals1"), gene_linear_plot2("Ralb"), 
          gene_linear_plot2("Dppa1"),  ncol = 4)
dev.off()


pdf(paste("GeneRidgePlot", Project, "16x_top_R_I_markers", "l2fc0.25", ".pdf", sep="_"), width=10, height=10)
par(bg=NA)
plot_grid(plt_ridges2("Phlda2"), plt_ridges2("Ncl"), plt_ridges2("Npm1"), plt_ridges2("Id2"), plt_ridges2("Paics"), 
          plt_ridges2("G3bp2"), plt_ridges2("Rpl4"), plt_ridges2("Cenpf"), plt_ridges2("Klf6"), 
          plt_ridges2("B2m"), plt_ridges2("Serpine2"), plt_ridges2("Krt18"), plt_ridges2("Acadl"), plt_ridges2("Lgals1"), 
          plt_ridges2("Ralb"), plt_ridges2("Dppa1"),  ncol = 4)
dev.off()


pdf(paste("GeneRidgePlot_w_linear", Project, "16x_top_R_I_markers", "l2fc0.25", ".pdf", sep="_"), width=10, height=10)
par(bg=NA)
plot_grid(plt_ridges_with_linear("Phlda2"), plt_ridges_with_linear("Ncl"), plt_ridges_with_linear("Npm1"), plt_ridges_with_linear("Id2"), 
          plt_ridges_with_linear("Paics"),  plt_ridges_with_linear("G3bp2"), plt_ridges_with_linear("Rpl4"), plt_ridges_with_linear("Cenpf"), 
          plt_ridges_with_linear("Klf6"),   plt_ridges_with_linear("B2m"), plt_ridges_with_linear("Serpine2"), plt_ridges_with_linear("Krt18"),
          plt_ridges_with_linear("Acadl"), plt_ridges_with_linear("Lgals1"),  plt_ridges_with_linear("Ralb"), plt_ridges_with_linear("Dppa1"),  ncol = 4)
dev.off()










