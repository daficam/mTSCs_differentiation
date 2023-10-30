#!/storage/Software/packages/anaconda3/bin/Rscript

rm(list=ls())
gc()
options(bitmapType='cairo')


library(cowplot)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(clustree)
library(harmony)
library(viridis)
library(colorspace)
library(optrees)
library(deldir)
library(alphahull)
library(RColorBrewer)

options(future.globals.maxSize = 4000 * 1024^2)

baseDir <- "/storage/CTR-Projects/CTR_DropSeq/NewMegatron/NewAlignments_50bp"
setwd(baseDir)
Project <- "COLOUR_ASSIGNMENT"






message("+-------------------------------------------------------------------------------")
message("+                        Plotting UMAP - custom                                 ")
message("+-------------------------------------------------------------------------------")

### load matrix.umap

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

reso                <- "SCT_snn_res.1.0"
clust               <- dropseq.integrated@meta.data %>% dplyr::select((SCT_snn_res.1))  # 

colnames(clust)     <- c("Cluster")
matrix.umap$Cluster <- clust$Cluster
matrix.umap <- matrix.umap[!grepl("LNA", matrix.umap$Experiment),] 

head(matrix.umap,2)
unique(matrix.umap$Cluster)
matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22","23", "24", "25")) # "23", "24", "25", "26" ))
matrix.umap$Cluster <- factor(matrix.umap$Cluster, levels = c("0","1","2","3","4","5","6","7","8", "9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22","23", "24", "25","26","27")) # "23", "24", "25", "26" ))


message("----------------------Cluster annotation -----------------------")

matrix.umap$Cluster <- as.numeric(as.character(matrix.umap$Cluster))
Clusters <- unique(matrix.umap$Cluster)
#Clusters <- as.numeric(levels(matrix.umap$Cluster))
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


message("----------------------umap.Cluster-----------------------")



vtess <- deldir::deldir(centers[,2], centers[,3], rw = c(range(dat[,1]), range(dat[,2])))
#plot(matrix.umap[,1], matrix.umap[,2], col = "gray", pch = 16, xlab = "UMAP 1", ylab = "UMAP 2")
#plot(vtess, wlines="triang", wpoints="dummy", number=FALSE, add=TRUE, lty=1, col = "black")
#points(centers[,2], centers[,3], cex = 3, pch = 21, bg = "white")
#text(centers[,2], centers[,3], labels = centers$Clusters, font = 2)

dat.tmp = centers
dist.tmp = dist(dat.tmp)

tri <- vtess$delsgs
tri$dist = tri$ind1
for(i in 1:nrow(tri)) 
{   
  tri$dist[i] = as.matrix(dist.tmp)[tri$ind1[i], tri$ind2[i]]
}
arcs = tri[,5:7]
colnames(arcs) = c("source", "target", "weight")
arcs = as.matrix(arcs)

BF_dists = list()
dists = vector()
for(k in 0:max(unique(Clusters)) ) {
  BF_dists[[k+1]] = optrees::getShortestPathTree(0:max(unique(Clusters)), arcs, source.node = k, algorithm = "Bellman-Ford", directed = F, show.data = F, show.graph = F, show.distances = F)
  dists[k+1] = sum(BF_dists[[k+1]]$distances) 
}

initial = which(dists == min(dists))-1

set.seed(82)



message("--------------------choose colour palette-----------------------")

nice_cols <-  c("#4d004b","#810f7c", "#88419d", "#bfd3e6","#8c96c6","#0868ac","#4eb3d3","#d9f0a3","#f768a1","#fa9fb5","#fde0dd","#d4b9da","#238443","#78c679","#004529","#081d58")
clcols <- nice_cols
pie(rep(1,length(nice_cols)), col=(nice_cols))

rainbow_cols <- rainbow(20)
n <- 20
pie(rep(1,n), col=(rainbow_cols))
#clcols <- rainbow_cols[c(1,3,4,5,6,8,12,14,15,16,17,20,21,22,23)]
clcols <- append(clcols, rainbow_cols[c(1,10,12,16,18,3,4,7)])
length(clcols)

n <- 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
pie(rep(1,n), col=(col_vector[1:50]))

clcols <- append(clcols, col_vector[c(4,7,8)])
length(clcols)
pie(rep(1,26), col=(clcols))

#clcols <- c("#141031","#99d992","#13f91d","#3107db","#d9bfdd","#300fd3","#19bb5f","#a30ef4","#d72ea0","#99a6aa","#5f3c58","#93f667","#2580a8","#b9ddfb","#7570ab","#08a022","#2e5ea0","#5ea9f6","#b96fa5","#35220d","#6b7b92","#f1ea64","#6ac573","#b447c1","#56434b","#d7eff9")

clcols <- clcols[-c(2,13,7,5,8,11)]

nice_cols <-  c("#4d004b","#810f7c", "#88419d", "#bfd3e6","#8c96c6","#0868ac","#4eb3d3","#d9f0a3","#f768a1","#fa9fb5","#fde0dd","#d4b9da","#238443","#78c679","#004529","#081d58")
clcols <- nice_cols

rainbow_cols <- rainbow(20)
n <- 20
pie(rep(1,n), col=(rainbow_cols))
#clcols <- rainbow_cols[c(1,3,4,5,6,8,12,14,15,16,17,20,21,22,23)]
clcols <- append(clcols, rainbow_cols[c(1,10,12,16,18,3,4)])
length(clcols)

n <- 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
pie(rep(1,n), col=(col_vector[1:50]))

clcols <- append(clcols, col_vector[c(4,7,8)])
length(clcols)
pie(rep(1,26), col=(clcols))

clcols <- append(clcols, c("#333399", "#cc99ff"))
pie(rep(1,28), col=(clcols))
length(clcols)

###
pie(rep(1,length(clcols)), col=(clcols))
clcols <- clcols[-c(3,4,7,10,13)]
pie(rep(1,length(clcols)), col=(clcols))
clcols <- clcols[-8] #21!!
####

clcols <- clcols[c(1:length(levels(dat$Cluster)))]


colordf = data.frame(levels(dat$Cluster), rep("gray", length(levels(dat$Cluster))), stringsAsFactors =  F)
colnames(colordf) = c("cluster", "color")
rownames(colordf) <- colordf$cluster
colordf[initial, 2] = sample(clcols, 1)
plot(matrix.umap[,1], matrix.umap[,2], col = colordf[matrix.umap[,"Cluster"],2], pch = 16, xlab = "UMAP 1", ylab = "UMAP 2")

nextcols = as.numeric(names(as.matrix(dist.tmp)[initial,][order(as.matrix(dist.tmp)[initial,])][2:ncol(as.matrix(dist.tmp))]))
nextcols

dist_to_initial = vector()
for(i in 1:length(clcols)) dist_to_initial[i] <- colorscience::deltaE2000(grDevices::convertColor(colorspace::hex2RGB(colordf[initial, 2])@coords, from = "sRGB", to = "Lab"),  grDevices::convertColor(colorspace::hex2RGB(clcols[(i)])@coords, from = "sRGB", to = "Lab") )


names(dist_to_initial) = clcols
dist_to_initial = dist_to_initial[order(dist_to_initial, decreasing = T)]
dist_to_initial = dist_to_initial[1:(length(dist_to_initial)-1)]
colordf[(nextcols+1),2] = names(dist_to_initial)

#plot(dat[,1], dat[,2], col = colordf[dat[,"Cluster"],2], pch = 16, xlab = "UMAP 1", ylab = "UMAP 2")

# Check how the assignment performed using Voronoi diagrams:
#plot.new()
#vtiles = tile.list(vtess)
#plot(vtiles, fill = colordf$color)
#points(dat[,1], dat[,2], col = "black", bg =  colordf[dat[,"Cluster"],2], pch = 21, xlab = NA, ylab = NA, cex = 0.7)
#points(centers[,"coord_x"], centers[,"coord_y"], cex = 3, pch = 21, bg = "white")
#text(centers[,"coord_x"], centers[,"coord_y"], labels = centers[,"Clusters"], font = 2)



library(spatstat)

cps <- ppp(centers[,"coord_x"], centers[,"coord_y"],window=owin(range(dat[,1]), range(dat[,2])))

sharededge <- function(X) {
  verifyclass(X, "ppp")
  Y <- X[as.rectangle(X)]
  dX <- deldir(Y)
  DS <- dX$dirsgs
  xyxy <- DS[,1:4]
  names(xyxy) <- c("x0","y0","x1","y1")
  sX <- as.psp(xyxy,window=dX$rw)
  marks(sX) <- 1:nobjects(sX)
  sX <- sX[as.owin(X)]
  tX <- tapply(lengths.psp(sX), marks(sX), sum)
  jj <- as.integer(names(tX))
  ans <- data.frame(ind1=DS[jj,5], 
                    ind2=DS[jj,6], 
                    leng=as.numeric(tX))
  return(ans)
}

shared_edge_lengths <- sharededge(cps)

# We now check which colors are within each shared edge, and their ΔE distances:
edgecolors = shared_edge_lengths
edgecolors$col1 = colordf[edgecolors$ind1,2]
edgecolors$col2 = colordf[edgecolors$ind2,2]
edgecolors$ind1 <- edgecolors$ind1-1
edgecolors$ind2 <- edgecolors$ind2-1
edgecolors$DeltaE.dist = apply(edgecolors, 1, function (x) colorscience::deltaE2000(grDevices::convertColor(colorspace::hex2RGB(x[4])@coords, from = "sRGB", to = "Lab"), grDevices::convertColor(colorspace::hex2RGB(x[5])@coords, from = "sRGB", to = "Lab")))

#By plotting all ΔE distances we see which two adjacent tessels have the lowest ΔE. For this palette and with this combination the distance values are not always optimal (defining “optimal” as above the arbitrary threshold of 20), but are always fairly distinguishable:

#plot(1:nrow(edgecolors), rep(0, nrow(edgecolors)), cex = 0, ylim = c(-1, max(edgecolors$DeltaE.dist)+10), xaxt = "n", xlab = NA, ylab = "DeltaE color distance")
#segments(x0 = 1:nrow(edgecolors), y0 = rep(0, nrow(edgecolors)), x1 = 1:nrow(edgecolors), y1 =edgecolors$DeltaE.dist)
#points(1:nrow(edgecolors), rep(0, nrow(edgecolors)), cex = 4, pch = 15, col = edgecolors$col1)
#points(1:nrow(edgecolors), edgecolors$DeltaE.dist, cex = 4, pch = 15, col = edgecolors$col2)
#axis(1, at = 1:nrow(edgecolors), labels = paste(edgecolors$ind1, edgecolors$ind2, sep = "-"), las = 2)


# 6: Cycle through all the colors to find the optimal starting point

colorassignments = list()
listcolordf = list()
for(j in 1:length(clcols)){   
  colordf = data.frame(0:(length(clcols)-1), rep("gray", length(clcols)), stringsAsFactors =  F)
  colnames(colordf) = c("cluster", "color")
  rownames(colordf) <- colordf$cluster
  colordf[ initial, 2] = clcols[j]
  dist_to_initial = vector()
  for(i in 1:length(clcols)) 
  {
    dist_to_initial[i] = colorscience::deltaE2000(grDevices::convertColor(colorspace::hex2RGB(colordf[initial, 2])@coords, from = "sRGB", to = "Lab"), grDevices::convertColor(colorspace::hex2RGB(clcols[i])@coords, from = "sRGB", to = "Lab"))
  }
  names(dist_to_initial) = clcols
  dist_to_initial = dist_to_initial[order(dist_to_initial, decreasing = T)]
  dist_to_initial = dist_to_initial[1:(length(dist_to_initial)-1)]
  #colordf[nextcols,2] = names(dist_to_initial)
  colordf[(nextcols+1),2] = names(dist_to_initial)
  edgecolors = shared_edge_lengths
  edgecolors$col1 = colordf[edgecolors$ind1,2]
  edgecolors$col2 = colordf[edgecolors$ind2,2]
  edgecolors$DeltaE.dist = apply(edgecolors, 1, 
                                 function (x) colorscience::deltaE2000(grDevices::convertColor(colorspace::hex2RGB(x[4])@coords, from = "sRGB", to = "Lab"), grDevices::convertColor(colorspace::hex2RGB(x[5])@coords, from = "sRGB", to = "Lab")))
  colorassignments[[j]] = edgecolors
  listcolordf[[j]] = colordf
}
mindists = unlist(lapply(colorassignments, function (x) min(x$DeltaE.dist[x$DeltaE.dist > 0])))

best.assignment = which(mindists == max(mindists))
if(length(best.assignment) > 1) best.assignment = best.assignment[1]

best.colordf = listcolordf[[best.assignment]]

plot(dat[,1], dat[,2], col = best.colordf[dat[,"Cluster"],2], pch = 16, xlab = "UMAP 1", ylab = "UMAP 2")
points(centers[,"coord_x"], centers[,"coord_y"], cex = 3, pch = 21, bg = "white")
text(centers[,"coord_x"], centers[,"coord_y"], labels = centers[,"Clusters"], font = 2)

best.colordf$color






















