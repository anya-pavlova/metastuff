#################################
###       beta-diversity      ###
#################################

ff <- sp_all[c(tags_vil, tags_control),]
ff <- ff[, which(colMaxs(ff) > 0)]

#rownames(ff) <- paste(meta_ourn[rownames(ff), "id_timed"],  meta_ourn[rownames(ff), "FIO"], sep=" |") 

# BC dist
dd <- bcdist(ff)

# UniFrac
#w_unifrac.m <- read.table("data/weighted_unifrac_merged_otu_tables.txt", header=T, row.names=1, sep="\t", colClasses = "character", as.is = T)
#colnames(w_unifrac.m) <- str_replace(colnames(w_unifrac.m), pattern = "^X", replacement = "")
#w_unifrac.m <- w_unifrac.m[c(tags_ea1, tags_ea2, tags_ab1, tags_ab2), c(tags_ea1, tags_ea2, tags_ab1, tags_ab2)]
#dd <- as.dist(w_unifrac.m)

mdd <- data.matrix(dd)
hr <- hclust(dd, method="ward")
groups <- cutree(hr, k=2)

if(proj_pref == "EE_KA") {
 print("CLUSTER SWAPPING - for illustration!!!")
 tmp <- groups
 tmp[which(groups == 1)] <- 2
 tmp[which(groups == 2)] <- 1
 groups <- tmp 
}
  

