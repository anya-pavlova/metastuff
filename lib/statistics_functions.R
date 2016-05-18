StatResult <- function(feat, group1, group2, maxpv, pairedt, outdir, nameplot, test){  
  if (test == TRUE){
    w <- wilcox_feats_both_dirs(feat, group1, group2, maxpv, pairedt)}
  else {w <- ttest_feats_both_dirs(feat, group1, group2, maxpv, pairedt)}
  
    w.greater <- as.data.table(w$greater, keep.rownames = T)
    setnames(w.greater,c("bacteria","pv"))
    w.greater <- w.greater[pv<=0.9]
    setkey(w.greater,'pv')
    w.greater$pv <- round(w.greater$pv, digits=5)
    w.greater$mediana <- round(w$greater_median_FC, digits=5)
    
    w.less <- as.data.table(w$less, keep.rownames = T)
    setnames(w.less,c("bacteria","pv"))
    w.less <- w.less[pv<=0.9]
    setkey(w.less,'pv')
    w.less$pv <- round(w.less$pv, digits=5)
    w.less$mediana <- round(w$less_median_FC, digits=5)
  
  outdir= '/home/anna/metagenome/16s_study/out_ee_ka/stat/'
  write.table(w.greater, 
              file = paste0(outdir, nameplot, "_stat_greatere.csv"), 
              sep = "\t", row.names = F, quote = F, col.names = T)
  write.table(w.less, 
              file = paste0(outdir, nameplot, "_stat_less.csv"), 
              sep = "\t", row.names = F, quote = F, col.names = T)
  
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
#T-test

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

