#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()

# run it on R.v4!!! :
# module load R/4.0.2

library(mrtree)
library(Seurat)
set.seed(42)



#dataset <- "HARMONY.orig_500g"
#dataset <- "HARMONY_I_500g"
dataset <- "HARMONY_R_500g"

if (dataset == "HARMONY.orig_500g"){
  Project <- "CTR_dscj1_HARMONY_orig_500g"
  reso <- "SCT_snn_res.0.8"
  baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY.orig_400g_500g/HARMONY.orig_500g"
  setwd(baseDir)
  matrix.su <- readRDS("HARMONY_matrix.su_umap_n.neigh25_repuls2_spread2L_harmony_ORIG_500g.Rds")
  clusters <- c(0:20) # for together dataset! - cluster 20 may need to be removed !
} else if (dataset == "HARMONY_I_500g"){
  Project <- "CTR_dscj1_HARMONY_I_500g"
  reso <- "SCT_snn_res.0.8"
  baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I"
  setwd(baseDir)
  matrix.su <- readRDS("harmony_500g_matrix.su_I.Rds")
  clusters <- c(0:16) 
  clust_order <-c("9",  "16", "11", "5", "1",  "15", "10", "12", "2",  "3",  "6",  "7",  "13", "8",  "14", "4" ,"0") # for I 500g res0.8
} else if (dataset == "HARMONY_R_500g"){
  Project <- "CTR_dscj1_HARMONY_R_500g"
  reso <- "SCT_snn_res.0.8"
  baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp/HARMONY_SEPARATE_R_I"
  setwd(baseDir)
  matrix.su <- readRDS("harmony_500g_matrix.su_R.Rds")
  clusters <- c(0:15) 
  clust_order <-c("14", "2",  "15", "7",  "13", "6",  "3",  "0",  "12", "5",  "10", "1",  "4",  "11", "9",  "8") # for R 500g res0.8
} else { print("please specify correct dataset!!")}



# create fer labels:::
Counts <- GetAssayData(matrix.su, assay= "SCT", slot= "data")
varfeats <- VariableFeatures(matrix.su)
Counts <- Counts[varfeats,]
Counts <- as.matrix(Counts)
Metadata <- matrix.su@meta.data
ref.labels <- Metadata$SCT_snn_res.0.8
head(ref.labels) # "Astrocytes" "Astrocytes" "Astrocytes" "Astrocytes" "Astrocytes"

# specify the resolution parameters
# resolutions = seq(0.1, sqrt(3), 0.05)^2
# alternatively and preferrably, we provide a sampling tool to sample resolution parameters to uniformly cover different scales
A = seurat_get_nn_graph(counts=Counts, metadata=Metadata, npc=10)
resolutions = modularity_event_sampling(A=A, n.res=30, gamma.min=0.01, gamma.max=2.5) # sample based on the similarity matrix
# clustering using Suerat 
seurat.out = sc_clustering.seurat(counts=Counts, resolutions=resolutions, metadata=Metadata, npcs=10,
                                  min.cells=0, min.features=0, scale.factor=10000, return.seurat.object=TRUE,
                                  vars.to.regress=NULL, find.variable.features=FALSE, verbose=FALSE)

# initial cluster tree from Seurat flat clustering
pdf(paste(Project, dataset, "seurat_flat_clustering", "1.pdf", sep="_"), onefile=FALSE, width=10, height=10) 
par(bg=NA)
plot_clustree(labelmat=seurat.out$seurat.clusters, prefix ='SCT_snn_res.', ref.labels = ref.labels, plot.ref = F)
dev.off()


# Then we apply MRtree to ubtain the hierarchical cluster tree, visualized using a dendrogram, with a pie chart on each tree node detailing the cluster composition given the known true labels.

out = mrtree(seurat.out$obj, consensus=F, augment.path=F, n.cores = 10)


pdf(paste(Project, dataset, "mrtree_hierarchicalClustree", "1.pdf", sep="_"), onefile=FALSE, width=10, height=10) 
par(bg=NA)
plot_tree(labelmat=out$labelmat.mrtree, ref.labels=ref.labels, plot.piechart = T,
          node.size = 0.4, tip.label.dist = 10, bottom.margin=30 )
dev.off()




# We evaluate te per-resolution clustering performance with a novel index adapted from Adjusted Rand Index to accrount for te bias for resolution.

ks.flat = apply(out$labelmat.flat, 2, FUN=function(x) length(unique(x)))
ks.mrtree = apply(out$labelmat.mrtree, 2, FUN=function(x) length(unique(x)))
amri.flat = sapply(1:ncol(out$labelmat.flat), function(i) AMRI(out$labelmat.flat[,i], ref.labels)$amri)
amri.flat = aggregate(amri.flat, by=list(k=ks.flat), FUN=mean)
amri.recon = sapply(1:ncol(out$labelmat.mrtree), function(i) AMRI(out$labelmat.mrtree[,i], ref.labels)$amri)
df = rbind(data.frame(k=amri.flat$k, amri=amri.flat$x, method='Seurat flat'), 
           data.frame(k=ks.mrtree, amri=amri.recon, method='MRtree'))
ggplot2::ggplot(data=df, aes(x=k, y=amri, color=method)) + geom_line() + theme_bw()

# We calcuate the similarity between the initial flat clustering and MRtree clusters across scales. Lower similarity indicates the selected clustering algorithm is not able to generate stabel clusters at the specific resolution. In this case stability drops steeply when $k>8$.

stab.out = stability_plot(out)
stab.out$plot




























