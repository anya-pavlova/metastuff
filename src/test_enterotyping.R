############################################################
#
# Cluster analysis of metagenomes (aka enterotyping).
# Code adapted from http://enterotype.embl.de/enterotypes.html
# Test on the 16S data from the original publication on enterotypes.
#
############################################################

#rm(list=ls(all=TRUE)); gc()
#setwd("~/metagenome/_external/biotyper/")

require(ade4)
require(clusterSim)

# set some functions

# Jensen-Shannon dissimilarity metric
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
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
  print(percent)
  return(Matrix_1)
}

# load test data (n=154 metagenomes)
# just a relative abundance table
load("data/test_data/Qin_el_al_enterotyping_data//Titanium16S.RData")

# put it into variable
data <- Titanium16S

# compute pairwise dissimilarity
data.dist=dist.JSD(data)


#nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

# find the optimal number of clusters (3 - for Qin et al data)
nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

# the cluster number that has maximum CH index is the winner (3 - here)
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

# just set it to 3
data.cluster=pam.clustering(data.dist, k=3)

# here are the clustering results - which sample belongs to which cluster
# this is what is interesting to us!
# I.e^ what is your cluster BEFORE, do you switch AFTER, how many people (in %) switch cluster like you, etc.
cbind(rownames(as.matrix(data.dist)), data.cluster)


# get silhouette width for the observed clustering (?)
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])

# denoise
data.denoized=noise.removal(data, percent=0.01)

# perform the between-class analysis
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 

# plot the result 
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)

