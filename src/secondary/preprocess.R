##################################################
# read OTU taxonomic assignment (Greengenes DB) ##
##################################################

taxotu <- read.table("data/97_otu_taxonomy.txt", header=F, row.names=1, sep="\t", colClasses = "character")

# TODO: read the global register of all metagenomes
#meta <- read.table("data/.txt", header=T, row.names=NULL, sep="\t", colClasses = "character")

#задаем выбор данных по нужным нам названиям
tags_vil <- sort(rownames(sp)[grep("^V", rownames(sp))])

# crop data from run
otu <- otu[tags_vil,]
otup <- otup[tags_vil,]
g <- g[tags_vil,]
fam <- fam[tags_vil,]
sp <- sp[tags_vil,]


# TODO - for base report - reads quality stats (results of Borya quality script )

# how many identified reads per sample
nreads <- rowSums(otu)
sort(nreads)
plot(sort(nreads))
MIN_NUM_READS <- 35000
good <- names(which(nreads >= MIN_NUM_READS))
# crop good covered
i <- which(tags_vil %in% good)
tags_vil <- tags_vil[i]

# crop data from run
otu <- otu[tags_vil,]
otup <- otup[tags_vil,]
g <- g[tags_vil,]
fam <- fam[tags_vil,]
sp <- sp[tags_vil,]


# compare ctrl and case coverage
boxplot(rowSums(otu_ctrl), rowSums(otu[tags_vil,]))

# examine Enterobacteriaceae content
#t <- fam[tags_vil, grep("Enterobacteriaceae", colnames(fam))]    
#t[,order(colMaxs(t))]


fam_all <- unite_feat_matrices(fam,fam_ctrl)
g_all <- unite_feat_matrices(g,g_ctrl)
sp_all <- unite_feat_matrices(sp,sp_ctrl)
otu_all <- unite_feat_matrices(otu,otu_ctrl)
otup_all <- unite_feat_matrices(otup,otup_ctrl)

tags_ctrl <- sort(rownames(sp_ctrl))

ff <- sp_all[c(tags_vil, tags_ctrl),]
ff <- ff[, which(colMaxs(ff) > 0)]


