slopegraph <- function(feats, feat_name, tags1, tags2, starts_red, my_ylab = "% of total bacterial abundance      ") #my_filename, my_width, my_height, 
{ 
#  feats <- al
 # feat_name <- "Alpha-diversity"
  #tags1 <- c(tags_ee1[1:length(tags_ee1)-1], tags_ka1)
  #tags2 <- c(tags_ee2[1:length(tags_ee2)-1], tags_ka2)
  #starts_red <-tags_ka1
  #my_ylab=""
  
  #tags_ee1 <- tags1
  #tags_ee2 <- tags2
  months<-24
  my_col <- rep("blue", length(tags1))
  names(my_col) <- tags1
  if(length(starts_red) > 0)
  {
    my_col[intersect(starts_red, tags1)] <- "red"
  }
  #print(my_col)
  year1 <- feats[tags1, feat_name] #gc(1338229205,5212325386,31725112511)
  year3 <- feats[tags2, feat_name] #c(1372425378,8836570075,49574919628)
  group <- tags1 #c("Group C", "Group B", "Group A")
  #print("1")
  a<-data.frame(year1,year3,group,my_col)
  #print("1.5")
  #print(a)
  #l11<-paste(a$group,comma_format()(round(a$year1)),sep="\n")#paste(a$group,comma_format()(round(a$year1/(3600*24*30.5))),sep="\n")
  l11<-round(a$year1,2)
 # print("1.51")
  #l13<-paste(a$group,comma_format()(round(a$year3)),sep="\n")
  l13<-round(a$year3,2)
 # print("1.52")
  #l13<-rep("", length(a$year3)) #round(a$year3,2)
  p<-ggplot(a) + geom_segment(data = a, aes(x=rep(0, length(year1)),xend=rep(24,length(year1)),
                                  y=year1,yend=year3),size=0.75,colour=my_col)
  
  ymax <- max(as.vector(feats[c(tags1, tags2), feat_name]))
  
  #print(identical(feat_name, colnames(feats)[grep("dolichum", colnames(feats))]))             
  #test <- max(as.vector(feats[c(tags_ab1, tags_ab2), feat_name]))
  #test <- max(as.vector(feats[c(tags_ab1, tags_ab2), grep("dolichum", colnames(feats))]))
  #test <- max(as.vector(ff[c(tags_ab1, tags_ab2), grep("dolichum", colnames(ff))]))  
  #print(ymax)
  #print(test)
 
 
  p<-p + theme(panel.background = element_blank())
  p<-p + theme(panel.grid=element_blank())
  p<-p + theme(axis.ticks=element_blank())
  p<-p + theme(axis.text=element_blank())
  p<-p + theme(panel.border=element_blank())
  p<-p + xlab("") + ylab(my_ylab) #+ +opts(axis.title.y = theme_text(vjust=0.6))
  p<-p + theme(axis.title.y=element_text(vjust=-3))
  p<-p + xlim((0-12),(months+12))
  p<-p + ylim(0,ymax)    #ylim(0,(1.2*(max(a$year3,a$year1))))
  
 # print("1.54")
  p<-p + geom_text(label=l13, y=a$year3, x=rep.int(months,length(a$year3)),hjust=-0.2,size=3.5)
 # print("1.541")
  p<-p + geom_text(label=l11, y=a$year1, x=rep.int( 0,length(a$year3)),hjust=1.2,size=3.5)
 # print("1.542")
  
 #p<-p + geom_text(label="Before", x=0,    y=(1.1*(max(a$year3,a$year1))),hjust= 1.2,size=5)
 p<-p + geom_text(label="", x=0,    y=6,hjust= 1.2, size=5)
 
 # print("1.543")
  #p<-p + geom_text(label="After", x=months,y=(1.1*(max(a$year3,a$year1))),hjust=-0.1,size=5)
 p<-p + geom_text(label="", x=months,y=6,hjust=-0.1,size=5)
 
 # print("1.544")
  p<-p + ggtitle(feat_name) + theme(legend.position="none")
 # print("1.545")
  #ggsave(filename=my_filename, p, device="pdf", width=my_width, height=my_height)
  return(p) 
}  




slopegraph_alpha <- function(feats, feat_name, tags1, tags2, starts_red, my_ylab = "% of total bacterial abundance      ") #my_filename, my_width, my_height, 
{ 
  #  feats <- al
  # feat_name <- "Alpha-diversity"
  #tags1 <- c(tags_ee1[1:length(tags_ee1)-1], tags_ka1)
  #tags2 <- c(tags_ee2[1:length(tags_ee2)-1], tags_ka2)
  #starts_red <-tags_ka1
  #my_ylab=""
  
  #tags_ee1 <- tags1
  #tags_ee2 <- tags2
  months<-24
  my_col <- rep("blue", length(tags1))
  names(my_col) <- tags1
  if(length(starts_red) > 0)
  {
    my_col[intersect(starts_red, tags1)] <- "red"
  }
  #print(my_col)
  year1 <- feats[tags1, feat_name] #gc(1338229205,5212325386,31725112511)
  year3 <- feats[tags2, feat_name] #c(1372425378,8836570075,49574919628)
  group <- tags1 #c("Group C", "Group B", "Group A")
  #print("1")
  a<-data.frame(year1,year3,group,my_col)
  #print("1.5")
  #print(a)
  #l11<-paste(a$group,comma_format()(round(a$year1)),sep="\n")#paste(a$group,comma_format()(round(a$year1/(3600*24*30.5))),sep="\n")
  l11<-round(a$year1,2)
  # print("1.51")
  #l13<-paste(a$group,comma_format()(round(a$year3)),sep="\n")
  l13<-round(a$year3,2)
  # print("1.52")
  #l13<-rep("", length(a$year3)) #round(a$year3,2)
  p<-ggplot(a) + geom_segment(data = a, aes(x=rep(0, length(year1)),xend=rep(24,length(year1)),
                                            y=year1,yend=year3),size=1.75,colour=my_col)
  
  #plot(p)
  
  
  # print("1.53")
  p<-p + theme(panel.background = element_blank())
  p<-p + theme(panel.grid=element_blank())
  p<-p + theme(axis.ticks=element_blank())
  p<-p + theme(axis.text=element_blank())
  p<-p + theme(panel.border=element_blank())
  p<-p + xlab("") + ylab(my_ylab) #+ +opts(axis.title.y = theme_text(vjust=0.6))
  p<-p + theme(axis.title.y=element_text(vjust=-3))
  p<-p + xlim((0-12),(months+12))
  p<-p + ylim(0,80)    #ylim(0,(1.2*(max(a$year3,a$year1))))
  
  # print("1.54")
  p<-p + geom_text(label=l13, y=a$year3, x=rep.int(months,length(a$year3)),hjust=-0.2,size=3.5)
  # print("1.541")
  p<-p + geom_text(label=l11, y=a$year1, x=rep.int( 0,length(a$year3)),hjust=1.2,size=3.5)
  # print("1.542")
  
  #p<-p + geom_text(label="Before", x=0,    y=(1.1*(max(a$year3,a$year1))),hjust= 1.2,size=5)
  p<-p + geom_text(label="", x=0,    y=6,hjust= 1.2, size=5)
  
  # print("1.543")
  #p<-p + geom_text(label="After", x=months,y=(1.1*(max(a$year3,a$year1))),hjust=-0.1,size=5)
  p<-p + geom_text(label="", x=months,y=6,hjust=-0.1,size=5)
  
  # print("1.544")
  p<-p + ggtitle(feat_name) + theme(legend.position="none")
  # print("1.545")
  #ggsave(filename=my_filename, p, device="pdf", width=my_width, height=my_height)
  return(p)
  
  #dev.off()
  #grid.arrange(p, ncol=1)
  #print(p)
  #print(p)
  #print("1.56")
}  

read_qiime_sum_feats <- function(filename)
{
  f <- t(read.table(filename, header=T, row.names=1, sep="\t"))
  rownames(f) <- gsub("X", "", rownames(f))
  f <- data.matrix(f)
  f <- 100 * f / rowSums(f)
  f
}

read_qiime_otu_table_no_tax <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t")
  # drop taxonomy
  f <- f[, -ncol(f)]
  f <- t(f)
  rownames(f) <- gsub("X", "", rownames(f))
  f
}

read_qiime_single_alpha_rar <- function(filename)
{   
  f <- read.table(filename, header=T, row.names=1, sep="\t", colClasses = "character")
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



#функция для склейки матриц разного размера
unite_feat_matrices <- function(t1, t2)
{
  uc <- sort(union(colnames(t1), colnames(t2)))
  ur <- c(rownames(t1), rownames(t2))
  t <- matrix(0, nrow = nrow(t1) + nrow(t2), ncol = length(uc))
  colnames(t) <- uc
  rownames(t) <- ur
  t[rownames(t1), colnames(t1)] <- t1
  t[rownames(t2), colnames(t2)] <- t2
  identical(t1, t[rownames(t1), colnames(t1)])
  identical(t2, t[rownames(t2), colnames(t2)])
  t
}

produce_barplots <- function(feat, tags, min_perc, separ, my_filename, my_width, my_height, rand_seed = 108, ...)
{  
  set.seed(rand_seed)
  f1 <- feat[tags,]
  f1 <- f1[,which(apply(f1, 2, max) >= min_perc)]
  colnames(f1) <- paste(separ, unlist(data.frame(strsplit(colnames(f1), separ))[2,]), sep="")
  rownames(f1) <- paste(meta_ourn[rownames(f1), "id_timed"],  meta_ourn[rownames(f1), "UID"], sep=" |") 
  df1 <- melt(f1)
  names(df1) <- c("sample", "Taxon", "abundance")
  print(getwd())
  print(my_filename)
  cairo_pdf(my_filename, width=my_width, height=my_height) #, ...)
  p <- ggplot(df1,aes(x=sample,y=abundance, fill=Taxon))+
    geom_histogram(stat="identity",colour = "black")+
    scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
    #theme(legend.text = element_text(colour="blue", size = 12))+
    ylab("Относит. представленность, %")+xlab("Образцы")+
    theme_bw()+ 
    theme(axis.text.x = element_text(angle = 90))    
    #theme(legend.position=element_textvjust())
  print(p)
  dev.off()
}

produce_barplots_ylimed <- function(feat, tags, min_perc, separ, my_filename, my_width, my_height, rand_seed = 108, my_ylim, ...)
{  
  set.seed(rand_seed)
  f1 <- feat[tags,]
  f1 <- f1[,which(apply(f1, 2, max) >= min_perc)]
  colnames(f1) <- paste(separ, unlist(data.frame(strsplit(colnames(f1), separ))[2,]), sep="")
  rownames(f1) <- paste(meta_ourn[rownames(f1), "id_timed"],  meta_ourn[rownames(f1), "UID"], sep=" |") 
  df1 <- melt(f1)
  names(df1) <- c("sample", "Taxon", "abundance")
  print(getwd())
  print(my_filename)
  cairo_pdf(my_filename, width=my_width, height=my_height) #, ...)
  p <- ggplot(df1,aes(x=sample,y=abundance, fill=Taxon))+
    geom_histogram(stat="identity",colour = "black")+
    scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
    #theme(legend.text = element_text(colour="blue", size = 12))+
    ylab("Относит. представленность, %")+xlab("Образцы")+
    theme_bw()+ 
    theme(axis.text.x = element_text(angle = 90))    +
    scale_y_continuous(limits = my_ylim)
  #theme(legend.position=element_textvjust())
  print(p)
  dev.off()
}

produce_barplots_pooled_diff <- function(feat, tags1, tags2, featnames1, featnames2, featnames3, featnames4, my_filename_inc, my_filename_dec, my_width, my_height, rand_seed = 108, leg_font_size = 6, ...)
{ 
  #     rand_seed <- 108
   #    feat <- foo
    #   tags1 <- tags_ea1
     #  tags2 <- tags_ea2
    #   featnames1 <- c("f__Lachnospiraceae; g__Lachnospira", "f__Lachnospiraceae; g__Roseburia", "f__; g__")
     #  featnames2 <- c("f__Lachnospiraceae; g__Lachnospira", "f__Lachnospiraceae; g__Roseburia")
      #  featnames3 <- c("f__Lachnospiraceae; g__Lachnospira", "f__Lachnospiraceae; g__Roseburia", "f__; g__") 
      #  featnames4 <- c("f__Lachnospiraceae; g__Lachnospira", "f__Lachnospiraceae; g__Roseburia")
      # leg_font_size <- 6
      
  set.seed(rand_seed)
  
  # reducing taxa
  if(length(featnames1) > 0)
  {
    f1 <- feat[c(tags1, tags2), featnames1, drop=F]    
    if(ncol(f1) > 1) {
      s1 <- colSums(f1[tags1, ])/length(tags1)
      s2 <- colSums(f1[tags2, ])/length(tags2)
      s2[setdiff(featnames1, featnames2)] <- 0
      df1 <- melt(cbind(s1, s2))
      names(df1) <- c("Taxon", "time", "abundance")
    } else {
      s1 <- as.vector(sum(f1[tags1, ])/length(tags1))
      s2 <- as.vector(sum(f1[tags2, ])/length(tags2))
      names(s1) <- featnames1
      names(s2) <- featnames2
      s2[setdiff(featnames1, featnames2)] <- 0
      df1 <- melt(cbind(s1, s2))        
      names(df1) <- c("Taxon", "time", "abundance")
    }   
    
    print(my_filename_dec)
    file.remove(my_filename_dec, showWarnings=F)
    cairo_pdf(my_filename_dec, width=my_width, height=my_height)
    pa <- ggplot(df1, aes(x=time, y=abundance, fill=Taxon))+
      geom_histogram(stat="identity",colour = "black")+    
      scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
      ylab("Средняя относит. представленность по группе, %")+xlab("Временная точка")+
      theme_bw()+ 
      theme(axis.text.x = element_text(angle = 90))   +
      scale_x_discrete(labels=c("до курса", "после курса")) + 
      theme(#legend.position="bottom",
            #legend.direction="vertical",
            #legend.key.size=unit(0.1, "cm"),
            legend.text=element_text(size=leg_font_size))
    print(pa)
    dev.off()
  }    
  
  # increasing taxa
  if(length(featnames3) > 0)
  {
    f1 <- feat[c(tags1, tags2), featnames3, drop=F]
    if(ncol(f1) > 1) {
      s1 <- colSums(f1[tags1, ])/length(tags1)
      s2 <- colSums(f1[tags2, ])/length(tags2)
      s2[setdiff(featnames3, featnames4)] <- 0
      df1 <- melt(cbind(s1, s2))        
      names(df1) <- c("Taxon", "time", "abundance")
    } else {
      s1 <- as.vector(sum(f1[tags1, ])/length(tags1))
      s2 <- as.vector(sum(f1[tags2, ])/length(tags2))
      names(s1) <- featnames2
      names(s2) <- featnames2
      s2[setdiff(featnames3, featnames4)] <- 0
      df1 <- melt(cbind(s1, s2))        
      names(df1) <- c("Taxon", "time", "abundance")
    }    
    print(my_filename_inc)
    file.remove(my_filename_inc, showWarnings=F)
    cairo_pdf(my_filename_inc, width=my_width, height=my_height)
    pb <- ggplot(df1, aes(x=time, y=abundance, fill=Taxon))+
      geom_histogram(stat="identity",colour = "black")+    
      scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
      ylab("Средняя относит. представленность по группе, %")+xlab("Временная точка")+
      theme_bw()+ 
      theme(axis.text.x = element_text(angle = 90))   +
      scale_x_discrete(labels=c("до курса", "после курса")) +
      theme(#legend.position="bottom",
        #legend.direction="vertical",
        #legend.key.size=unit(0.1, "cm"),
        legend.text=element_text(size=leg_font_size))
    print(pb)
    dev.off()
  }
}









produce_barplots_pooled_diff_yet_another <- function(feat, tags1, tags2, featnames_inc, featnames_dec, my_filename_inc, my_filename_dec, my_width, my_height, rand_seed = 108, leg_font_size = 6, ...)
{ 
#     rand_seed <- 108
#     feat <- foo
#     tags1 <- tags_ea1
#     tags2 <- tags_ea2
#     featnames_inc <- c()
#     featnames_dec <- c("f__; g__", "f__Verrucomicrobiaceae; g__Akkermansia")
#                        
#     my_filename_inc <- paste("graphs/", proj_pref, "/diff_multibarplots_ea_inc_g.pdf", sep="") 
#     my_filename_dec <- paste("graphs/", proj_pref, "/diff_multibarplots_ea_dec_g.pdf", sep="")
#     my_width <- 6
#     my_height <- 8
#    leg_font_size <- 6
#   
  
  
  
  set.seed(rand_seed)  
  # reducing taxa
  if(length(featnames_dec) > 0)
  {
    f1 <- feat[c(tags1, tags2), featnames_dec, drop=F]    
    if(ncol(f1) > 1) {
      s1 <- colSums(f1[tags1, ])/length(tags1)
      s2 <- colSums(f1[tags2, ])/length(tags2)
      df1 <- melt(cbind(s1, s2))
      names(df1) <- c("Taxon", "time", "abundance")
    } else {
      s1 <- as.vector(sum(f1[tags1, ])/length(tags1))
      s2 <- as.vector(sum(f1[tags2, ])/length(tags2))      
      names(s1) <- featnames_dec
      names(s2) <- featnames_dec
      df1 <- melt(cbind(s1, s2))
      names(df1) <- c("Taxon", "time", "abundance")
    }   
    print(my_filename_dec)
    file.remove(my_filename_dec, showWarnings=F)
    cairo_pdf(my_filename_dec, width=my_width, height=my_height)     # , order = -as.numeric(Taxon))
    pa <- ggplot(df1, aes(x=time, y=abundance, fill=Taxon, order = -as.numeric(Taxon)))+
      geom_histogram(stat="identity", colour = "black")+    
      scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
      ylab("Средняя относит. представленность по группе, %")+xlab("Временная точка")+
      theme_bw()+ 
      theme(axis.text.x = element_text(angle = 90))   +
      scale_x_discrete(labels=c("до курса", "после курса")) + 
      theme(legend.text=element_text(size=leg_font_size))
    print(pa)
    dev.off()
  }    
  
  # increasing taxa
  if(length(featnames_inc) > 0)
  {
    f1 <- feat[c(tags1, tags2), featnames_inc, drop=F]
    if(ncol(f1) > 1) {
      s1 <- colSums(f1[tags1, ])/length(tags1)
      s2 <- colSums(f1[tags2, ])/length(tags2)
      df1 <- melt(cbind(s1, s2))        
      names(df1) <- c("Taxon", "time", "abundance")
    } else {
      s1 <- as.vector(sum(f1[tags1, ])/length(tags1))
      s2 <- as.vector(sum(f1[tags2, ])/length(tags2))
      names(s1) <- featnames_inc
      names(s2) <- featnames_inc      
      df1 <- melt(cbind(s1, s2))        
      names(df1) <- c("Taxon", "time", "abundance")
    }    
    print(my_filename_inc)
    file.remove(my_filename_inc, showWarnings=F)
    cairo_pdf(my_filename_inc, width=my_width, height=my_height)
    pb <- ggplot(df1, aes(x=time, y=abundance, fill=Taxon, order = -as.numeric(Taxon)))+
      geom_histogram(stat="identity",colour = "black")+    
      scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
      ylab("Средняя относит. представленность по группе, %")+xlab("Временная точка")+
      theme_bw()+ 
      theme(axis.text.x = element_text(angle = 90))   +
      scale_x_discrete(labels=c("до курса", "после курса")) +
      theme(#legend.position="bottom",
        #legend.direction="vertical",
        #legend.key.size=unit(0.1, "cm"),
        legend.text=element_text(size=leg_font_size))
    print(pb)
    dev.off()
  }
  
  
  
}




produce_barplots_simple <- function(in_filename, out_filename, my_width, my_height, rand_seed = 108, leg_font_size = 6, ...)
{
  set.seed(rand_seed)
  f <- read.table(in_filename, header=F, row.names=1, sep="\t")
  f <- f[which(!is.na(f[,1])),,drop=F]
  f <- f[which(f[,1] != 0),,drop=F]
  a <- rep(1.0, nrow(f))
  b <- 1/f[,1]
  #b[which(f[,1] == "Inf"),] <- 0
  f2 <- as.matrix(cbind(a, b))
  rownames(f2) <- rownames(f)
  colnames(f2) <- c("t1", "t2")
  f <- f2
  df1 <- melt(f)
  names(df1) <- c("Taxon", "time", "abundance")
  file.remove(out_filename, showWarnings=F)
  cairo_pdf(out_filename, width=my_width, height=my_height)
  pa <- ggplot(df1, aes(x=time, y=abundance, fill=Taxon, order = -as.numeric(Taxon)))+
    geom_histogram(stat="identity",colour = "black")+    
    scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
    ylab("Изменение, раз")+xlab("Временная точка")+
    theme_bw()+ 
    theme(axis.text.x = element_text(angle = 90))   +  
    theme(legend.text=element_text(size=leg_font_size))
  print(pa)
  dev.off()
}
  
  
  


convert_taxa_diff_to_table <- function(t_d)
{
  a <- data.frame(names(t_d$greater), as.vector(t_d$greater), round(t_d$greater_median_FC, 2), stringsAsFactors = F)
  b <- data.frame(names(t_d$less), as.vector(t_d$less), round(t_d$less_median_FC, 2), stringsAsFactors = F)
  a <- a[order(a[,3], decreasing=T),]
  colnames(a) <- c("taxon", "FDR adj. p-value", "fold change")
  b <- b[order(b[,3], decreasing=F),]  
  colnames(b) <- c("taxon", "FDR adj. p-value", "fold change")
  tt_d <- rbind(a,b)                  
  tt_d
}

convert_taxa_diff_to_table_DIV <- function(t_d)
{
  t_d <- w  
  
  a <- data.frame(names(t_d$greater), as.vector(t_d$greater), round(t_d$greater_median_FC, 2), stringsAsFactors = F)
  b <- data.frame(names(t_d$less), as.vector(t_d$less), round(t_d$less_median_FC, 2), stringsAsFactors = F)
  a <- a[order(a[,3], decreasing=T),]
  colnames(a) <- c("taxon", "FDR adj. p-value", "fold change")
  b <- b[order(b[,3], decreasing=F),]  
  colnames(b) <- c("taxon", "FDR adj. p-value", "fold change")
  div <- data.frame('-', 1, 1)
  colnames(div) <- c("taxon", "FDR adj. p-value", "fold change")
  tt_d <- rbind(a, div, b)
  tt_d
#     
#   tt_d <- rbind.data.frame(a, c("-", 1, 1), div)
# 
#   div <- data.frame(`taxon` = '-', 
#            `FDR adj. p-value` = 0,
#            `fold change` = 0)
#   names(div) <- names(a)
#   q <- rbind(a, div, b)
# 
#  q <- rbind(a, div, b)
#   q[,2]  
# 
#   div[,2]
#   a[,2]
#   
#   tt_d[,2]
}





delta_scatter_graph <- function(feats, feat_name)
{
  b <- feat_name #colnames(ff)[1]
  delta_ea <- ff[tags_ea2, b] - ff[tags_ea1, b]
  delta_ab <- ff[tags_ab2, b] - ff[tags_ab1, b]
  temp <- rbind(data.frame(class=rep(name_ea, length(delta_ea)),delta=delta_ea), data.frame(class=rep(name_ab, length(delta_ab)),delta=delta_ab))
  #temp$alp <- rep(1, nrow(temp))
  #temp$alp[which(abs(temp$delta) < 0.5)] <- 0    
  ggplot(temp, aes(x=class, y=delta)) + geom_jitter(position=position_jitter(0.05), alpha=0.4) + labs(title=b, x="Group", y = "Delta %") + theme_classic()   +  theme(plot.title = element_text(size=7))
  
}


produce_barplots_bottom_legend <- function(feat, tags, min_perc, separ, my_filename, my_width, my_height, rand_seed = 108, ncol=5, ltextsize=10, lsize=0.3, ...)
{  
  set.seed(rand_seed)
  f1 <- feat[tags,]
  f1 <- f1[,which(apply(f1, 2, max) >= min_perc)]
  colnames(f1) <- paste(separ, unlist(data.frame(strsplit(colnames(f1), separ))[2,]), sep="")
  rownames(f1) <- paste(meta_ourn[rownames(f1), "id_timed"],  meta_ourn[rownames(f1), "UID"], sep=" |") 
  df1 <- melt(f1)
  names(df1) <- c("sample", "Taxon", "abundance")
  print(getwd())
  print(my_filename)
  cairo_pdf(my_filename, width=my_width, height=my_height) #, ...)
  p <- ggplot(df1,aes(x=sample,y=abundance, fill=Taxon))+
    geom_histogram(stat="identity",colour = "black")+
    scale_fill_manual(values=sample(rainbow(n=nrow(df1))))+
    #theme(legend.text = element_text(colour="blue", size = 12))+
    ylab("Относит. представленность, %")+xlab("Образцы")+
    theme_bw()+ 
    theme(axis.text.x = element_text(angle = 90))+
    theme(legend.position="bottom",
          legend.direction="vertical",
          legend.key.size= unit(lsize, "cm"),
          legend.text = element_text(size=ltextsize))+
    guides(fill=guide_legend(ncol=ncol))
    #theme(legend.position=element_textvjust())
  print(p)
  dev.off()
}

