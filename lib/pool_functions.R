GetTopInfo <- function(group, col.prefix, top.limit=3){
  # Get top N bacteria from group
  #
  # group - matrix. Format:
  #          _bac1_|_..._|_bacK_
  # sample1 |      |     |
  # ...     |      |     |
  # sampleM |      |     |
  #
  # col.prefix - string, column prefix
  #
  # result - 
  # 1) data.table with columns:
  # | <col.prefix>sample | <col.prefix>1 | ... | <col.prefix>N | <col.prefix>1p | ... | <col.prefix>Np |
  #
  # columns <col.prefix>N - bacteria names
  # columns <col.prefix>Np - bacteria abundance percent
  # 2) list: top bacterias by position

  if(top.limit == -1){ top.limit <-  length(colnames(group))}
  top.by.sample <- alply(.data = group, 
                         .margins = 1, 
                         .fun = function(x) sort(x, decreasing = TRUE))
  names(top.by.sample) <- rownames(group)
  top.by.sample.table <- as.data.table(ldply(.data = top.by.sample, 
                                             .fun = function(x) c(names(x[1:top.limit]), x[1:top.limit])))
  setnames(top.by.sample.table, 
           c(paste(col.prefix, "sample", sep=""), 
             paste(col.prefix, 1:top.limit, sep=""), 
             paste(col.prefix, 1:top.limit, "p", sep="")))
  
  top.pool.by.position <- alply(.data=1:top.limit, 
                                .margins=1, 
                                .fun=function(i) unique(rapply(top.by.sample, 
                                                          function(x) names(x[i]))))
  pool.top <- as.data.table(unlist(lapply(top.by.sample, function(x) names(x)[1:top.limit])))
  setnames(pool.top, c("bacteria"))
  pool.top <- pool.top[,.N,by=bacteria]
  setorder(pool.top, -N)
   
  list(top.by.sample.table=top.by.sample.table, 
       top.pool.by.position=top.pool.by.position,
       top.by.sample=top.by.sample,
       pool.top=pool.top)
}

GetPool <- function(matrix, thresh.percent){
  threshold.percent <- thresh.percent
  is.in.pool <- apply(matrix, 2, FUN = function(x) sum(x > threshold.percent))
  pool <- as.data.table(is.in.pool, keep.rownames = T)
  setnames(pool, c("bacteria", "count"))
  setorder(pool, -count)
  pool[count > 0]
}

GetTopPlot <- function(sample.top){
  top.before <- as.data.table(sample.top, keep.rownames = T)
  setnames(top.before, c("bacteria", "percent"))
  top.before$bacteria <- make.short.name(top.before$bacteria)
  ggplot(top.before[1:20], aes(x = reorder(bacteria, percent), y = percent)) + 
    geom_bar(stat="identity", width = 0.5) + 
    theme_bw() + 
    xlab('bacteria') + 
    ylab('proportion') + 
    theme(text = element_text(size=5)) +
    coord_flip()
}

GroupPool <- function(featurse, meta.data.g1, meta.data.g2, group.name, thresh.percent){
  pool.before <- GetPool(featurse[rownames(featurse) %in% meta.data.g1,], thresh.percent)
  pool.after <- GetPool(featurse[rownames(featurse) %in% meta.data.g2,], thresh.percent)
  
  write.table(pool.before, 
              file = paste0(pool_directory,group.name, "_species_pool_before.csv"), 
              sep = "\t", row.names = F, quote = F, col.names = T)
  
  write.table(pool.after, 
              file = paste0(pool_directory, group.name,"_species_pool_after.csv"), 
              sep = "\t", row.names = F, quote = F, col.names = T)
}
