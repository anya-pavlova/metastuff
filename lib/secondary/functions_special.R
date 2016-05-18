# sum by group (GO, KEGG Orthology, COG)
sumPOGByGroup <- function(tt, anno) 
{  
  q <- as.data.frame(tt)
  q$sample <- rownames(tt)
  q.m <- melt(q, id.vars = 'sample')
  q.mm <- merge(q.m, anno[, c("bgi_name", "group_name")], by.x = 'variable', by.y = 'bgi_name')
  #head(q.mm)  
  q.a <- aggregate(data = q.mm, value ~ group_name + sample, sum)
  #head(q.a)  
  q.c <- cast(q.a, sample ~ group_name, value = 'value')
  rownames(q.c) <- q.c[,"sample"]
  res <- q.c[,-1]  
  res  
}

# sum by group (GO, KEGG Orthology, COG) and binarize
sumPOGByGroupBin <- function(tt, anno) 
{  
  q <- as.data.frame(tt)
  q$sample <- rownames(tt)
  q.m <- melt(q, id.vars = 'sample')
  q.mm <- merge(q.m, anno[, c("bgi_name", "group_name")], by.x = 'variable', by.y = 'bgi_name')
  #head(q.mm)  
  q.a <- aggregate(data = q.mm, value ~ group_name + sample, sum)
  #head(q.a)  
  q.c <- cast(q.a, sample ~ group_name, value = 'value')
  rownames(q.c) <- q.c[,"sample"]
  res <- q.c[,-1]
  # binarize
  res <- ifelse(res == 0, 0, 1)
  res  
}

# find binary features rare in one set and frequent in another at the same time
specificBinFeats <- function(ff, b1, b2, thres_freq=0.8, thres_rare=0.1)
{
  # frequent in 1
  i1 <- which( colSums(ff[b1,]) >= thres_freq * length(b1))
  print(paste("num frequent in arg 1: ", length(i1), sep=" "))
  g1 <- colnames(ff)[i1]
  # rare in 2
  i2 <- which( colSums(ff[b2,]) <= thres_rare * length(b2))
  print(paste("num rare in arg 2: ", length(i2), sep=" "))
  g2 <- colnames(ff)[i2]
  # intersect
  choi <- intersect(g1, g2)
  print(paste("intersection: ", length(choi), sep=" "))
  choi
#   qq <- which(try3$sseqid %in% choi)
#   length(qq)
#   genes_evil_tab <- try3[qq, c("qseqid", "sseqid")]
#   write.table(genes_evil_tab, file="out/roots_of_evil.txt", quote=F, row.names=F, sep="\t")
#   sss <- sort(t(data.frame(strsplit(as.vector(genes_evil_tab[,1]), split="\\|"), stringsAsFactors=F)[4,]))
#   write.table(sss, file="out/roots_of_evil_add.txt", quote=F, row.names=F, sep="\t")
#   genes_evil <- genes_evil_tab$sseqid
}

# convert list of BGI genes IDs to our ortholgy groups
# WARNING! does not preserve order!
convertBGItoOGs_unordered <- function(list_bgi, panan)
{    
  ii <- which(panan$sseqid %in% list_bgi)
  sort(unique(panan[ii, "OG_joined"]))  
}






