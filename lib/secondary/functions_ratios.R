

# for paired metagenomes, get fold change values by dividing:
# importanly - first value by second
get_decrease_ratios <- function(feats, tags1, tags2)
{ 
  ff <- feats[c(tags1, tags2), ]  
  #colnames(ff) <- paste("g_", unlist(data.frame(strsplit(colnames(ff), "g_"))[2,]), sep="")
  ff <- ff[,which(colSums(ff) > 0)]
  ff <- ff + 1E-10 # pseudocount
  
  # ratio of increase
  #rs <- ff[tags2, ] / ff[tags1, ]  
  # OR   
  # ratio of decrease
  rs <- ff[tags1, ] / ff[tags2, ]
  
  rownames(rs) <- tags2
  rs
}



compare_decrease_ratios <- function(feat, tags_ea1, tags_ea2, tags_ab1, tags_ab2, cmp_prefix, pval_thresh, doAdjust, cut_point)
{  
  #feat <- g
  #pval_thresh <- 0.5
  #doAdjust <- T
  #cut_point <- "o__"  
  
  rs_ea <- get_decrease_ratios(ff, tags_ea1, tags_ea2)
  rs_ab <- get_decrease_ratios(ff, tags_ab1, tags_ab2)
  
  names_rs_com <- intersect(colnames(rs_ea), colnames(rs_ab))
  rs_com <- rbind(rs_ea[,names_rs_com], rs_ab[,names_rs_com])  
  
  difs <- c()
  for(tax in names_rs_com)
  {  
    # initialize flag 
    cur_dif <- 1
    
    #print(tax)
    cmp_ea <- rs_com[tags_ea2, tax]; del <- which(is.na(cmp_ea)); if(length(del) > 0) { cmp_ea <- cmp_ea[-del]; }
    cmp_ab <- rs_com[tags_ab2, tax]; del <- which(is.na(cmp_ab)); if(length(del) > 0) { cmp_ab <- cmp_ab[-del]; }
    
    setdiff(tags_ea2, rownames(rs_com))
    
    #boxplot(cmp_ea, cmp_ab)
    
    tt <- c()
    
    #########################
    # for analysis of ratio of increase
    ##t <- try(t.test(cmp_ea, cmp_ab, alternative = "greater"), silent=TRUE)
    #t <- try(wilcox.test(cmp_ea, cmp_ab, alternative = "greater"), silent=TRUE)
    # OR
    # for analysis of ratio of decrease
    ##t <- try(t.test(cmp_ea, cmp_ab, alternative = "less"), silent=TRUE)
    t <- try(wilcox.test(cmp_ea, cmp_ab, alternative = "less"), silent=TRUE)
    ########################
    
    if (is(t, "try-error")) {    
      # print("ttest exception")    
      cur_dif <- -1
    }
    else {
      # print(t$p.value)
      cur_dif <- t$p.value
    }
    
#     # crop names
#     names(cur_dif) <- paste(cut_point, unlist(data.frame(strsplit(tax, cut_point))[2,]), sep="") 
#     # OR
#     # leave names
     names(cur_dif) <- tax

     difs <- append(difs, cur_dif)
  }
  difs <- data.frame(difs)
  colnames(difs) <- "p-value"
  
  # adjust p-values if flagged
  if(doAdjust)
  {
    p_adj <- p.adjust(difs[, "p-value"], method='fdr')
    difs[, "p-value"] <- p_adj
  }
    
  # filter by p-value
  difs <- difs[which(difs[, "p-value"] <= pval_thresh),,drop=F]

  if(nrow(difs) > 0)
  {
    # crop names
    rownames(difs) <- paste(cut_point, unlist(data.frame(strsplit(rownames(difs), cut_point))[2,]), sep="") 
  }

  print(difs)  
  write.table(difs, file=paste("out/", proj_pref, "/cmp_decrease_dif_", cmp_prefix, ".txt", sep=""), sep="\t", quote=F, row.names = T)  
  
}

