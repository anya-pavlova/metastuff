############################################################
#
# Cluster analysis of metagenomes (aka enterotyping).
# Code adapted from http://enterotype.embl.de/enterotypes.html
#
############################################################

library(ade4)
library(clusterSim)
source("lib/coocurance_function.R")

# palettes:
# library(RColorBrewer)
# display.brewer.all()
cluster.colors <- NA

# Jensen-Shannon dissimilarity metric
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

# partition around medoids clustering
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

# noise removal - drop low abundance taxa
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  return(Matrix_1)
}

GetClusters <- function(abundance.matrix, data.samples){
  data <- t(abundance.matrix)
  message('start clustering')
  # denoise
  data <- noise.removal(data, percent=0.01)
  
  # compute pairwise dissimilarity
  data.dist=dist.JSD(data)
  
  message('choosing K')
  # find the optimal number of clusters
  # k = 1 and k = 2 values give an error further
  nclusters=NULL
  for (k in 3:20) { 
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
  
  # the cluster number that has maximum CH index is the winner
  # plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
  k <- which.max(nclusters)
  message(paste('K is', k))
  
  data.cluster=pam.clustering(data.dist, k)
  
  # merge data.samples and cluster ids
  data.sample.cluster <- data.table(cbind(rownames(as.matrix(data.dist)), data.cluster))
  setnames(data.sample.cluster, c("sample","cluster"))
  cluster.before <- data.sample.cluster[sample %in% data.samples$before]
  setnames(cluster.before, c("sample", "cluster_before"))
  cluster.after <- data.sample.cluster[sample %in% data.samples$after]
  setnames(cluster.after, c("sample", "cluster_after"))

  data.samples <- merge(data.samples, cluster.before, by.x = "before", by.y = "sample")
  data.samples <- merge(data.samples, cluster.after, by.x = "after", by.y = "sample")
  data.samples$moved <- data.samples$cluster_before != data.samples$cluster_after
  
  message('between-class analysis')
  # perform the between-class analysis to plot clusters and find out drivers
  pca.result <- dudi.pca(data.frame(t(data)), scannf=F, nf=10)
  bca.result <- bca(pca.result, fac=as.factor(data.cluster), scannf=F, nf=k-1)
  
  # bca.result$tab contains info about drivers: which features contributed to cluster mostly
  # we take [1] - first driver
  drivers <- unlist(alply(bca.result$tab, 1, function(x) sort(x, decreasing = T)[1]))
  
  message('preparing clusters info')
  # now get cluster info
  cluster.info <- GetClusterInfo(data.samples)
  cluster.info$cluster.perc$driver <- names(drivers)
  
  message('finish clustering')
  return(list(data.sample.cluster.paired = data.samples, 
              data.sample.cluster = data.sample.cluster, 
              bca.result = bca.result,
              cluster.perc = cluster.info$cluster.perc,
              cluster.moves = cluster.info$movements,
              data.cluster = data.cluster))
}

PlotClusters <- function(bca.result, data.cluster, drivers){
  cluster.colors <<- brewer.pal(length(drivers), "Set1")
  # plot the result 
  s.class(bca.result$ls, 
          fac=as.factor(data.cluster), 
          grid=F, 
          label = drivers, 
          col = cluster.colors, 
          clabel = 2)
}

# plot before and after for one sample
PlotBeforeAfterClusters <- function(bca.result, data.sample.cluster, sample.before, sample.after){
  bca.coordinates.before <- bca.result$ls[rownames(bca.result$ls) == sample.before,]
  cluster.before <- as.integer(data.sample.cluster[sample == sample.before]$cluster)
  points(bca.coordinates.before$CS1, bca.coordinates.before$CS2, pch=24, col="black", bg=cluster.colors[cluster.before], cex=3)

  bca.coordinates.after <- bca.result$ls[rownames(bca.result$ls) == sample.after,]
  cluster.after <- as.integer(data.sample.cluster[sample == sample.after]$cluster)
  points(bca.coordinates.after$CS1, bca.coordinates.after$CS2, pch=22, col="black", bg=cluster.colors[cluster.after], cex=3)
}

GetClusterInfo <- function(data.sample.cluster.paired){
  # Gets info from cluster analysis

  # % of people for each cluster before and after
  cluster.perc.before <- data.sample.cluster.paired[,(.N/nrow(data.sample.cluster.paired))*100, by=cluster_before]
  cluster.perc.after <- data.sample.cluster.paired[,(.N/nrow(data.sample.cluster.paired))*100, by=cluster_after]
  cluster.perc <- merge(cluster.perc.before, cluster.perc.after, by.x='cluster_before', by.y='cluster_after')
  setnames(cluster.perc, c("cluster", "before", "after"))
  
  # % of people moved from which cluster to which cluster
  movements <- data.sample.cluster.paired[moved==T, 
                                          (.N/nrow(data.sample.cluster.paired))*100, 
                                          by=.(cluster_before, cluster_after)]
  setnames(movements, "V1", "perc")
  setorder(movements, -perc)
  
  return(list(cluster.perc=cluster.perc, movements=movements))
}