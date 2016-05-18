################################################
# functions                                    #
################################################

#loading libraries
library(futile.logger)
library(lint)
library(stringr)
library(ecodist)
library(data.table)
library(ape)
library(ggplot2)
library(gplots)
library(scales)
library(MASS)
library(matrixStats)
library(gridExtra)
library(grid)
library(phyloseq) #source("https://bioconductor.org/biocLite.R"), biocLite("phyloseq")
library(reshape)
library(asbio)
library(vegan)

##################
# Functions: 
# UniteMatrix, ReadIni, WriteTable, QIIME functions,     

# UniteMatrices 
# input: two matrix  
# output: united matrix
# TODO - uncertain - sensitive to arguments order - see !!!
UniteMatrices <- function(t1, t2)
{
  if (all(is.na(t1))){return(t2)}
  if (all(is.na(t2))){return(t1)}
  uc <- sort(union(colnames(t1), colnames(t2)))
  ur <- c(rownames(t1), rownames(t2))
  t <- matrix(0, nrow = nrow(t1) + nrow(t2), ncol = length(uc))
  colnames(t) <- uc
  rownames(t) <- ur
  # !!!
  t[rownames(t1), colnames(t1)] <- t1
  t[rownames(t2), colnames(t2)] <- t2
  identical(t1, t[rownames(t1), colnames(t1)])
  identical(t2, t[rownames(t2), colnames(t2)])
  t
}

# ReadIni
# input: .ini file 
# output: list
ReadIni <- function(IniFilename) 
{ 
  connection <- file(IniFilename) 
  Lines  <- readLines(connection) 
  close(connection) 
  
  Lines <- chartr("[]", "==", Lines)  # change section headers 
  
  connection <- textConnection(Lines) 
  d <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE) 
  close(connection) 
  
  L <- d$V1 == ""                    # location of section breaks 
  d <- subset(transform(d, V3 = V2[which(L)[cumsum(L)]])[1:3], 
              V1 != "") 
  
  ToParse  <- paste("INI.list$", d$V3, "$",  d$V1, " <- '", 
                    d$V2, "'", sep="") 
  
  INI.list <- list() 
  eval(parse(text=ToParse)) 
  
  return(INI.list) 
} 

# WriteTable
# input: some structure
# output: txt file with data frome input file
WriteTable <- function (someStructure, outdir, type)
{
  filename<-paste(outdir, "/", type, '.csv', sep="")
  write.table(someStructure, filename, quote=T, sep='\t')
}

# QIIME functions 
# input: tables frome qiime
# output: matrix with data from input file
ReadQiimeSumFeats <- function(filename)
{
  f <- t(read.table(filename, header=T, row.names=1, sep="\t"))
  rownames(f) <- gsub("X", "", rownames(f))
  f <- data.matrix(f)
  f <- 100 * f / rowSums(f)
  f
}

ReadQiimeOtuTableNoTax <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t")
  # drop taxonomy
  f <- f[, -ncol(f)]
  f <- t(f)
  rownames(f) <- gsub("X", "", rownames(f))
  f
}

ReadQiimeSingleAlphaRar <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t", 
                  colClasses = "character")
  f <- f[,-c(1,2)]
  colnames(f) <- gsub("X", "", colnames(f))
  f <- t(f)
  colnames(f) <- "
  chao1"
  f <- f[which(f[,"chao1"] != "n/a"),,drop=F]
  f <- data.frame(f, stringsAsFactors = F)
  f[,"chao1"] <- as.numeric(f[,"chao1"])
  f
}

# loading alpha diversity
LoadAlphaDiv <- function(inpdir)
{
  inpdir
  AlphaDivTbl <- (read.table (inpdir, row.names=1,  header=T,sep="\t",
                              stringsAsFactors=F ))
  AlphaDivTbl <- t(AlphaDivTbl)
  AlphaDivTbl <- as.data.frame(AlphaDivTbl[-c(1, 2),])
  rownames(AlphaDivTbl)<-gsub("X","", rownames(AlphaDivTbl))
  colnames(AlphaDivTbl) <- c("AlphaDiversity")
  return(data.matrix(AlphaDivTbl))
}


# chooseTOPfeature 
# choose top features with total % of abundance across all samples
# input: feature vectors (table), percantage
# output: table with chosen features
chooseTOPfeature<-function(data,perc)
{
  summof<-sum(colSums(data))
  sorted<-sort(colSums(data), decreasing = TRUE)
  count <-0
  num<-0
  name_list<-c()  
  for (i in sorted) #in sorted by total abundance features, take those making %
  {
    count <- count+((i/summof)*100)
    num <- num+1
    if (count > perc) 
    {
      break
    }
    else 
    {
      name_list <- append(name_list, names(sorted[num]))
    }    
  }  
  y <- data[,name_list]
  return(y)
}



# make  distance in Spearman correlation  - for internal use
# input: feature vectors (table)
# output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}

# good colours for heatplot
cols.gentleman <- function(ncol=500) {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, 'RdBu'))(ncol)
  return(rev(hmcol))
}


# make  distance in Spearman correlation  - for internal use
# input: feature vectors (table)
# output: class dist object
distSpear<-function (x) 
{
  result <- 1-cor(t(x), method='spearman', use='pair')
  as.dist(result)
}

# make JS distance - for internal use
# input: feature vectors (table)
# output: class dist object
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


# make all distances
# input: feature vectors (table), vector with names of metrics ('JS','Spear','Eu','Man','Can','BC') chosen
# output: a list with class dist objects, names of the objects - like metrics chosen
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


HeatMap <- function(features,outdir, plot_name )
{
  # features = ee.features$species.ee.1
  # outdir = '/home/anna/metagenome/16s_study/out_ee_ka/heatmap'
  # plot_name = 'TEST'
  ff_par_top <- features[, which(colMaxs(features) > 3)]
  ff_par_top <- ff_par_top[, order(colMaxs(ff_par_top), decreasing = F)]
  dd <- bcdist(ff_par_top)
  hr <- hclust(dd, method="ward.D")
  hc <- hclust(distSpear(t(ff_par_top)), method="ward.D")
  
  # final pre-print change names of heatmap rows and cols
  #colnames(ff_par_top) <- paste("o_", unlist(data.frame(strsplit(colnames(ff_par_top), "o_"))[2,]), sep="")
  # draw
  cairo_pdf(paste(outdir, '/HeatMap_', plot_name, '.pdf', sep=''), width=35, height=30)
  #cairo_pdf(paste("graphs/heat_fam.pdf", sep=""), width=15, height=10)
  heatmap.2(data.matrix(ff_par_top), 
            Rowv=as.dendrogram(hr),
            distfun = function(x) dist(x, method = "euclidean"),
            dendrogram="row",
            col=cols.gentleman(500), 
            #cexRow=1.2, cexCol=1.2, trace='none', scale="none", margin=c(15,20) )
            cexRow=1.2, cexCol=1.2, trace='none', scale="none", margin=c(25,20) )
  
  dev.off()
}  
################################################################################################### 
###################################################################################################


### Print Mean and Stdev for group of samples
PrintMeanAndStd<-function(fam.gen.spe.listTOP, min.trsh, group.name){
  
  MeansTOP <- lapply(fam.gen.spe.listTOP, colMeans) 
  SdTOP <- lapply(fam.gen.spe.listTOP, function(x){apply(x,2,sd)})
  
  for (i in 1:length(MeansTOP)){
    ResultMeanAndSd<-cbind(round(MeansTOP[[i]][which(MeansTOP[[i]]>=min.trsh)],digits=2), 
                           round(SdTOP[[i]][which(MeansTOP[[i]]>=min.trsh)], digits=2))
    ResultMeanAndSd <- as.data.table(ResultMeanAndSd, keep.rownames = T)
    colnames(ResultMeanAndSd) <- c('bacteria',"Mean, %", "Stdev, %")
    WriteTable (ResultMeanAndSd, "/home/anna/metagenome/16s_study/out_ee_ka", paste(group.name,"_tableTOP_",i)) 
  }
}

PrintMeanAndStd2<-function(fam.gen.spe.listTOP, min.trsh, pathway.out){
  MeansTOP <- lapply(fam.gen.spe.listTOP, colMeans) 
  SdTOP <- lapply(fam.gen.spe.listTOP, function(x){apply(x,2,sd)})
  for (i in 1:length(MeansTOP)){
    ResultMeanAndSd<-cbind(round(MeansTOP[[i]][which(MeansTOP[[i]]>=min.trsh)],digits=2), 
                           round(SdTOP[[i]][which(MeansTOP[[i]]>=min.trsh)], digits=2))
    colnames(ResultMeanAndSd)<-c("Mean, %", "Stdev, %")
    WriteTable (ResultMeanAndSd, pathway.out, paste("tableTOP_g2_",i)) 
  }
  ### TO DO: print for one sample
}

################################################################################
# Stat.tests functions                                                        ##
################################################################################
#ToDo: check that the distribution is normal?  
ttest_feats_both_dirs <- function(feat, group1, group2, maxpv=0.05, pairedt, ...)
{ 
  pg <- c() # p-values for greater
  pl <- c() # p-values for less
  fcg <- c() # median fold change for greater
  fcl <- c() # median fold change for less
  for(i in 1:ncol(feat))
  {  
    w <- t.test(feat[group1, i], feat[group2, i], paired=pairedt, alternative="greater", ...)    
    if(!is.na(w$p.value))
    {
      tmp <- w$p.value
      names(tmp) <- colnames(feat)[i]
      pg <- append(pg, tmp)    
      fcg <- append(fcg, median(feat[group1, i])/median(feat[group2, i]))
    }
    
    w <- t.test(feat[group1, i], feat[group2, i], paired=pairedt, alternative="less", ...)
    if(!is.na(w$p.value))
    {
      tmp <- w$p.value
      names(tmp) <- colnames(feat)[i]
      pl <- append(pl, tmp)
      fcl <- append(fcl, median(feat[group1, i])/median(feat[group2, i]))
    }
  }  
  pg_adj <- p.adjust(pg, method='fdr')
  pl_adj <- p.adjust(pl, method='fdr')      
  i <- which(pg_adj <= maxpv)
  pg_adj <- pg_adj[i]
  fcg <- fcg[i]
  i <- which(pl_adj <= maxpv)
  pl_adj <- pl_adj[i]  
  fcl <- fcl[i]
  res <- list(pg_adj, fcg, pl_adj, fcl)  
  names(res) <- c("greater", "greater_median_FC", "less", "less_median_FC")
  res
}

# test log fold change for each featurem then adjust p-values
# WARNING! I drop values where feat is zero in at least one of two samples
ttest_paired_log_fc <- function(feat, group1, group2, maxpv=0.05)
{   
  
  #feat=ff
  #group1=tags_sol_il[1,]
  #group2=tags_sol_il[2,]
  #i =1 
  
  pg <- c()
  pl <- c()  
  fcg <- c() # median fold change for greater
  fcl <- c() # median fold change for less  
  for(i in 1:ncol(feat))
  {
    #print(colnames(feat)[i])
    
    tt <- log(feat[group1, i] / feat[group2, i])
    # leave finite values
    tt <- tt[is.finite(tt)]
    # if one value is left, then not enough, skip it
    if(length(tt) <= 1)
    {
      next
    }
    w <- t.test(tt, mu=0, alternative="greater")    
    if(!is.na(w$p.value))
    {
      tmp <- w$p.value
      names(tmp) <- colnames(feat)[i]
      pg <- append(pg, tmp)    
      fcg <- append(fcg, median(tt))
    }    
    
    w <- t.test(tt, mu=0, alternative="less")
    if(!is.na(w$p.value))
    {
      tmp <- w$p.value
      names(tmp) <- colnames(feat)[i]
      pl <- append(pl, tmp)
      fcl <- append(fcl, median(tt))
    }
  }  
  #print("im out")
  pg_adj <- p.adjust(pg, method='fdr')
  pl_adj <- p.adjust(pl, method='fdr')      
  i <- which(pg_adj <= maxpv)
  pg_adj <- pg_adj[i]
  fcg <- fcg[i]
  i <- which(pl_adj <= maxpv)
  pl_adj <- pl_adj[i]  
  fcl <- fcl[i]
  res <- list(pg_adj, fcg, pl_adj, fcl)  
  names(res) <- c("greater", "greater_median_logFC", "less", "less_median_logFC")
  res
}

wilcox_feats_both_dirs <- function(feat, group1, group2, maxpv=0.05, pairedt)
{   
  pg <- c()
  pl <- c()  
  fcg <- c() # median fold change for greater
  fcl <- c() # median fold change for less
  for(i in 1:ncol(feat))
  {  
    w <- wilcox.test(feat[group1, i], feat[group2, i], paired=pairedt, alternative="greater")    
    if(!is.na(w$p.value))
    {
      tmp <- as.numeric(w$p.value)
      names(tmp) <- colnames(feat)[i]
      pg <- append(pg, tmp)    
      fcg <- append(fcg, as.numeric((median(feat[group1, i])/median(feat[group2, i]))))
    }
    
    w <- wilcox.test(feat[group1, i], feat[group2, i], paired=pairedt, alternative="less")
    if(!is.na(w$p.value))
    {
      tmp <- as.numeric(w$p.value)
      names(tmp) <- colnames(feat)[i]
      pl <- append(pl, tmp)
      fcl <- append(fcl, as.numeric(median(feat[group1, i])/median(feat[group2, i])))
    }
  }  
  pg_adj <- p.adjust(pg, method='fdr')
  pl_adj <- p.adjust(pl, method='fdr')      
  i <- which(pg_adj <= maxpv)
  pg_adj <- pg_adj[i]
  fcg <- fcg[i]
  i <- which(pl_adj <= maxpv)
  pl_adj <- pl_adj[i]  
  fcl <- fcl[i]
  res <- list(pg_adj, fcg, pl_adj, fcl)
  names(res) <- c("greater", "greater_median_FC", "less", "less_median_FC")
  res
}

# paired; output mean of deltas (not fold change of medians)
wilcox_feats_both_dirs_PAIRED_MEAN_DELTA <- function(feat, group1, group2, maxpv=0.05)
{   
  pairedt <- TRUE
  pg <- c()
  pl <- c()  
  fcg <- c() # median fold change for greater
  fcl <- c() # median fold change for less
  for(i in 1:ncol(feat))
  {  
    print(colnames(feat)[i])
    w <- wilcox.test(feat[group1, i], feat[group2, i], paired=pairedt, alternative="greater")    
    if(!is.na(w$p.value))
    {
      tmp <- as.numeric(w$p.value)
      names(tmp) <- colnames(feat)[i]
      pg <- append(pg, tmp)
      #fcg <- append(fcg, as.numeric((median(feat[group1, i])/median(feat[group2, i]))))
      tmp <- feat[group2, i] - feat[group1, i]
      print(tmp)
      fcg <- append(fcg, as.numeric(mean(tmp)))
    }
    
    w <- wilcox.test(feat[group1, i], feat[group2, i], paired=pairedt, alternative="less")
    if(!is.na(w$p.value))
    {
      tmp <- as.numeric(w$p.value)
      names(tmp) <- colnames(feat)[i]
      pl <- append(pl, tmp)
      #fcl <- append(fcl, as.numeric(median(feat[group1, i])/median(feat[group2, i])))
      tmp <- feat[group2, i] - feat[group1, i]
      print(tmp)
      fcl <- append(fcl, as.numeric(mean(tmp)))
    }
  }  
  pg_adj <- p.adjust(pg, method='fdr')
  pl_adj <- p.adjust(pl, method='fdr')      
  i <- which(pg_adj <= maxpv)
  pg_adj <- pg_adj[i]
  fcg <- fcg[i]
  i <- which(pl_adj <= maxpv)
  pl_adj <- pl_adj[i]  
  fcl <- fcl[i]
  res <- list(pg_adj, fcg, pl_adj, fcl)
  names(res) <- c("greater", "greater_mean_delta", "less", "less_mean_delta")
  res
}

wilcox_onesample_both_dirs <- function(feat, group, x, maxpv=0.05)
{   
  pg <- c() # p-values for greater
  pl <- c() # p-values for less
  fcg <- c() # median fold change for greater
  fcl <- c() # median fold change for less
  for(i in 1:ncol(feat))
  {  
    w <- wilcox.test(feat[group, i], mu = feat[x, i], alternative="greater")    
    if(!is.na(w$p.value))
    {
      tmp <- w$p.value
      names(tmp) <- colnames(feat)[i]
      pg <- append(pg, tmp)    
      fcg <- append(fcg, median(feat[group, i])/feat[x, i])
    }
    
    w <- wilcox.test(feat[group, i], mu = feat[x, i], alternative="less")
    if(!is.na(w$p.value))
    {
      tmp <- w$p.value
      names(tmp) <- colnames(feat)[i]
      pl <- append(pl, tmp)
      fcl <- append(fcl, median(feat[group, i])/feat[x, i])
    }
  }    
  pg_adj <- p.adjust(pg, method='fdr')
  pl_adj <- p.adjust(pl, method='fdr')      
  i <- which(pg_adj <= maxpv)
  pg_adj <- pg_adj[i]
  fcg <- fcg[i]
  i <- which(pl_adj <= maxpv)
  pl_adj <- pl_adj[i]  
  fcl <- fcl[i]
  res <- list(pg_adj, fcg, pl_adj, fcl)  
  names(res) <- c("greater", "greater_median_FC", "less", "less_median_FC")
  res
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
