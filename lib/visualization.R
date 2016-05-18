########################################################
## Visualization: MDS, HeatMaps, BoxPlot              ##
########################################################

flog.info("make MDS for species levels")

samp.id.all <-data.groups$sample_id
filtred.sample.id.all <- samp.id.all[samp.id.all %in% row.names(species)]# образцы из samp.id.all кт есть in species 
bact.for.samp.id.all <- species[filtred.sample.id.all, ]

dist.all<-bcdist(bact.for.samp.id.all)
dist.all.m <- as.matrix(dist.all)

bfr.samples.mtrx <- dist.all.m[data.samples$before, data.samples$before] 
dist.bfr.samples <- as.dist(bfr.samples.mtrx)
aftr.samples.mtrx <- dist.all.m[data.samples$after, data.samples$after] 
dist.aftr.samples <- as.dist(aftr.samples.mtrx)

setkey(data.groups, 'group_id')
meta.for.bfr.samples <- data.groups["BEFORE_TREAT"]
meta.for.aftr.samples <- data.groups["AFTER_TREAT"]
meta.for.existing.bact <- data.groups[data.groups$sample_id %in% filtred.sample.id.all,]

MDS(dist.all, pathway$outputdir$graphs.dir, meta.for.existing.bact, "group_id", 'species.tow.groups') 
MDS(dist.bfr.samples, pathway$outputdir$graphs.dir, meta.for.bfr.samples, "group_id", 'group.before')
MDS(dist.aftr.samples, pathway$outputdir$graphs.dir, meta.for.aftr.samples, "group_id", 'group.after' )

MDS<-function(dis, outdir, meta, colFact, plot_name, m_width = 7, m_height=5)
{ 
  myMDS <- isoMDS(dis, k=2)
  outdirMDS<-paste(outdir,'/MDS_',plot_name,'.pdf',sep='')
  outdirMDS<-file.path(outdirMDS)
  #without lables
  dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]), groups = meta[[colFact]], lab = ' ')
  #dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]),cols = meta[,colFact], lab = rownames(myMDS$points))
  plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=groups, label=lab)) + geom_point(size=3, alpha=0.7) +
            #geom_text(aes(x=X+2, y=Y, size=2)) + theme_bw() 
             theme_bw()+
             theme(legend.key = element_blank())
  ggsave(outdirMDS, plotMDS, width = m_width, height=m_height)  
}

IndividualMDS <-function(dis.all, pathway$outputdir$graphs.dir,
                         meta.for.existing.bact, "group_id", 'INDMDS', 
                         m_width = 7, m_height=5)

IndividualMDS <-function(dis, outdir, meta, colFact, plot_name, m_width = 7, m_height=5)
{
  dis=dist.all
  outdir= pathway$outputdir$graphs.dir
  meta= meta.for.existing.bact
  colFact="group_id"
  plot_name = "paired.samp.test"
  m_width = 7
  m_height=5
  
  
  myMDS <- isoMDS(dis, k=2)
  outdirMDS<-paste(outdir,'/MDS_',plot_name,'.pdf',sep='')
  outdirMDS<-file.path(outdirMDS)
  dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]), groups = meta[[colFact]], lab = rownames(myMDS$points))
  indiv.dfMDS <- dfMDS[which(dfMDS$lab %in% data.samples[which(data.samples[[1]] == '1.121.1'),]),]
  plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=groups, label=lab)) + geom_point(size=3,alpha=0.7) +  
    #geom_text(hjust=0, vjust=0, col="black")+
    #sample before
    geom_point(data = indiv.dfMDS[1,], size = 7, color='orchid4',alpha=0.9)+ 
    geom_text(data = indiv.dfMDS[1,], col="black", label="you before" )+
    #sample after
    geom_point(data = indiv.dfMDS[2,], size = 7, color='orchid4',alpha=0.9 )+ 
    geom_text(data = indiv.dfMDS[2,], col="black", label="you after" )+
    theme_bw()+
    theme(legend.key = element_blank())
  
  ggsave(outdirMDS, plotMDS, width = m_width, height=m_height)  
}


IndividualMDS <-function(dis, outdir, meta, colFact, plot_name, m_width = 7, m_height=5)
{
  
  myMDS <- isoMDS(dis, k=2)
  outdirMDS<-paste(outdir,'/MDS_',plot_name,'.pdf',sep='')
  outdirMDS<-file.path(outdirMDS)
  dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]),cols = meta[,colFact], lab = rownames(myMDS$points))
  indiv.dfMDS <-dfMDS[which(dfMDS$lab %in% totalTable$meta.pairedsamples[which(totalTable$meta.pairedsamples[,1] == "1.155.1"),]),]
  
  plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=factor(cols), label=lab)) + geom_point(size=3) + #geom_text(aes(x=X+2, y=Y, size=2)) + theme_bw() 
    geom_text(hjust=0, vjust=0, col="black")+theme_bw()+
    #sample before
    geom_point(data = indiv.dfMDS[1,], size = 12, color='darkorchid4')+ 
    geom_text(data = indiv.dfMDS[1,], col="black", label="you before" )+
    #sample after
    geom_point(data = indiv.dfMDS[2,], size = 12, color='darkorchid4')+ 
    geom_text(data = indiv.dfMDS[2,], col="black", label="you after" )  
  
  ggsave(outdirMDS, plotMDS, width = m_width, height=m_height)  
}

IndividualMDS(dist.all, pathway$group.one$outdir, meta.for.existing.bact, "quality_status", "MDS_individual")


#samp.id.all - all samples which is in meta table
samp.id.all <-totalTable$meta.samplesgroups[totalTable$meta.samplesgroups[,"group_id"] %in% c("BEFORE_TREAT","AFTER_TREAT"), 'sample_id' ]
#filtred.sample.id.all - sampls from samp.id.before which is in the totalTable$family
filtred.sample.id.all <- samp.id.all[samp.id.all %in% row.names(totalTable$family)]
#part from totalTable$family just for filtred.sample.id.all - samples
bact.for.samp.id.all <- totalTable$family[filtred.sample.id.all, ]
#meta tabel for samples which we have
meta.for.existing.bact <- totalTable$meta.samplesgroups[totalTable$meta.samplesgroups[, 'sample_id'] 
                                                        %in% row.names(bact.for.samp.id.all),]

flog.info("make heatmaps")
#Making HeatMap
HeatMap(ee.features$species.ee.1, 'out_ee_ka/heatmap/', 'family.ee.bfr')
HeatMap(feature.list.g1$family, pathway$outputdir$graphs.dir, 'family.avva.g1' )
HeatMap(feature.list.g1$genus, pathway$outputdir$graphs.dir, 'genus.avva.g1')
HeatMap(feature.list.g1$species, pathway$outputdir$graphs.dir, 'species.avva.g1')

HeatMap(feature.list.g2$family, pathway$outputdir$graphs.dir, 'family.avva.g2')
HeatMap(feature.list.g2$genus, pathway$outputdir$graphs.dir, 'genus.avva.g2')
HeatMap(feature.list.g2$species, pathway$outputdir$graphs.dir, 'species.avva.g2')

#Making HeatMap for TOP Features
HeatMap(feature.TOPlist.g1$family, pathway$outputdir$graphs.dir, 'family.avvaTOPfeatures.g1' )
HeatMap(feature.TOPlist.g1$genus, pathway$outputdir$graphs.dir, 'genus.avvaTOPfeatures.g1' )
HeatMap(feature.TOPlist.g1$species, pathway$outputdir$graphs.dir, 'species.avvaTOPfeatures.g1')

HeatMap(feature.TOPlist.g2$family, pathway$outputdir$graphs.dir, 'family.avvaTOPfeatures.g2' )
HeatMap(feature.TOPlist.g2$genus, pathway$outputdir$graphs.dir, 'genus.avvaTOPfeatures.g2' )
HeatMap(feature.TOPlist.g2$species, pathway$outputdir$graphs.dir, 'species.avvaTOPfeatures.g2')


#######################################################################
#MDS for family, genus, species
# for (i in names(totalTable)[1:3]){
#   samp.id.all <-totalTable$meta.samplesgroups[totalTable$meta.samplesgroups[,"group_id"] %in% c("BEFORE_TREAT","AFTER_TREAT"), 'sample_id' ]
#   filtred.sample.id.all <- samp.id.all[samp.id.all %in% row.names(totalTable[[i]])]# образцы из samp.id.before кт есть totalTable$family 
#   bact.for.samp.id.all <- totalTable[[i]][filtred.sample.id.all, ]
#   meta.for.existing.bact <- totalTable$meta.samplesgroups[totalTable$meta.samplesgroups[, 'sample_id'] %in% row.names(bact.for.samp.id.all),]
#   
#   samp.id.before <-totalTable$meta.samplesgroups[totalTable$meta.samplesgroups[,"group_id"] %in% "BEFORE_TREAT", 'sample_id' ]
#   filtred.sample.id.before <- samp.id.before[samp.id.before %in% row.names(totalTable[[i]])]# образцы из samp.id.before кт есть totalTable$family 
#   bact.for.samp.id.before <-totalTable[[i]][filtred.sample.id.before, ]
#   meta.for.existing.bact.bef <- totalTable$meta.samplesgroups[totalTable$meta.samplesgroups[, 'sample_id'] %in% row.names(bact.for.samp.id.before),]
#   
#   dist.all<-bcdist(bact.for.samp.id.all)
#   MDS(dist.all, pathway$outputdir$graphs.dir, meta.for.existing.bact,"group_id", i) 
#   
#   dist.bf<-bcdist(bact.for.samp.id.before)
#   MDS(dist.bf, pathway$outputdir$graphs.dir, meta.for.existing.bact.bef,"group_id", paste('befor_', i)) 
#   
# }
  
  

