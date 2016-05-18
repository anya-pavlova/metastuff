###make  distance in Spearman correlation  - for internal use
####input: feature vectors (table)
####output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}
#####edn


###make JS distance - for internal use
####input: feature vectors (table)
####output: class dist object
distJS<-function(x)
{
  result<-c()
  input<-c()
  covs<-x
  covs<-covs+0.000001 # add pseudocount
  covs_norm<-covs/rowSums(covs) #normalize input by 1
  input = t(covs_norm)
  heads<-list(colnames(input))
  result<-matrix(0,nrow = ncol(input), ncol = ncol(input),dimnames = c(heads, heads)) #make result distance table
  kld1<-0
  kld2<-0
  for (i in 1:ncol(input))
  {
    for (p in 1:ncol(input)) 
    {
      xj<-input[,i]
      yj<-input[,p]
      mj<-(xj+yj)/2
      kld1<-kld1 + xj*(log(xj/mj)) #compute Kullback-Leibler coefficients
      kld2<-kld2 + yj*(log(yj/mj))
      result[i,p]<-sqrt(sum((kld1/2) ,(kld2/2))) #compute Jenson-Shannon distance
      kld1<-0
      kld2<-0
    }  
  }
  js_d<-as.dist(result)
  return (js_d)
}
#######end

###make all distances
####input: feature vectors (table), vector with names of metrics ('JS','Spear','Eu','Man','Can','BC') chosen
####output: a list with class dist objects, names of the objects - like metrics chosen
allDist<-function(x, metric=c('JS','Spear','Eu','Man','Can','BC'))
{
  dists<-list()
  num<-0
  for (i in (metric)) 
  {
{
  num<-num+1
  if (i == 'JS')
  {
    dis<-distJS(x)
  }
  if (i == 'Spear')
  {
    dis<-distSpear(x)
  }
  if (i == 'Eu')
  {
    dis<-dist(x)
  }
  if (i == 'Man')
  {
    dis<-dist(x, method='manhattan')
  }
  if (i == 'Can')
  {
    dis<-dist(x, method='canberra')
  }
  if (i == 'BC')
  {
    dis<-bcdist(x)
  }
  dists[[num]]<-dis
}	
  }	
names(dists)<-metric
return (dists)
}
