#################################
###       beta-diversity      ###
#################################

ff <- sp_all[c(tags_vil, tags_ctrl),] #sp_all?
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
hr <- hclust(dd, method="ward.D") #The "ward" method has been renamed to "ward.D"; note new "ward.D2"
groups <- cutree(hr, k=2)

#makeMDS<-function(distObj, meta, colFact, symbFact, outdir)
makeMDS(dd, , , , /out)
