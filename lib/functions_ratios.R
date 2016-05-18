############################################################################
# functions for finding taxa that change in % between the two timepoints
# AT
############################################################################

# constants
MAX_RATIO <- 1e6 # if ratio > this value, round it to this to make ranked tests valid.


# For a pair of metagenomes, get fold change values, divide the changing into increasers and decreasers,
# leave only major taxa (with max. rel. abundance higher than threshold).
get.changed <- function(feats, tag1, tag2, min_max_perc=0.1, min_change_ratio=1.0)
{ 
  tf <- as.data.table(t(feats[c(tag1, tag2),]), keep.rownames = T)
  setnames(tf, c('rn', tag1, tag2), c('name', 'sample1', 'sample2'))
  tf[(sample2>0)&(sample1>0), ratio := sample1/sample2]
  tf[(sample2==0)&(sample1>0), ratio := MAX_RATIO]
  tf[(sample2>0)&(sample1==0), ratio := 1/MAX_RATIO]
  tf[(sample2==0)&(sample1==0), ratio := 0.0]
  
  tf.thresh <-tf[tf[,max(sample1, sample2) > min_max_perc, by = 1:nrow(tf)]$V1]
  
  # divide into increasers and decreasers
  decs <- tf.thresh[ratio > 1, c('ratio', 'name'), with=F]
  decs <- setorder(decs, 'ratio')
  decs <- decs[ratio > min_change_ratio]
  
  incs <- tf.thresh[(0 < ratio) & (ratio < 1), c('ratio', 'name'), with=F]
  incs[, ratio:=1/ratio]  
  incs <- setorder(incs, 'ratio')
  incs <- incs[ratio > min_change_ratio]
  
  # pack the results
  return(list("taxa_decreased"=decs, "taxa_increased"=incs))
}
