#EE and KA
rm(list=ls(all=TRUE)); gc()

current_directory <- "/home/anna/metagenome/16s_antib"
setwd(current_directory)

library(data.table)
library(reshape)
library(plyr)
library(RColorBrewer)
library(HMP)
library(grid)

source("lib/functions.R")
source("lib/read_functions.R")
source("lib/coocurance_function.R")
source("lib/pool_functions.R")
source("lib/rank_abundance_curve.R")
source("lib/picrust.R")
source("lib/cluster_analysis.R")
source("lib/functions_ratios.R") # analysis of differences in bacterial abundance
source("lib/diversity.R")
source("lib/choose_functions.R")
source("lib/load_functions.R")
source("functions_from_Vera/functions_power.R")
source("lib/statistics_functions.R")

# todo: bad
make.short.name <- function(vector.names) {
  vector.names.new <- str_replace(vector.names, ".*NN=", "")
  vector.names.new <- str_replace(vector.names.new, "\\|D=.*", "")
  vector.names.new <- str_replace(vector.names.new, ".*;", "")
  vector.names.new <- strsplit(vector.names.new, "_", fixed=T)
  vector.names.new <- sapply(vector.names.new, function(x) {paste(x[1], x[2])})
  return(vector.names.new)
}

path.to.config <- file.path(current_directory, 'src', 'pathway.ini')
pathway <- ReadIni(path.to.config) 
path.to.config <- file.path(current_directory, 'src', 'pathway_ee_ka.ini')
pathway.ee.ka <- ReadIni(path.to.config) 

# loading meta data for EE and KA groups
# LoadMeta
meta.list <- MetaLoading(pathway$meta.tables$meta.pairedsamples,
                         pathway$meta.tables$meta.samplesgroups)


#loading all data
name.family <- "otu_table_L4.txt"
name.genus <- "otu_table_L5.txt"
name.species <- "otu_table_L6.txt"

totalTable <- LoadProjects(projects.directory.path=pathway$main$project.folder, 
                           file.name.family = name.family,
                           file.name.genus = name.genus,
                           file.name.species = name.species,
                           folder = 'summarize_taxa_perc')

#meta data
data.project.samples <- meta.list$data.samples.ee.ka.pair
data.samples.ee.ka.pair <- meta.list$data.samples.ee.ka.pair

#choosing project samples
family <- ChooseProjectSamples(totalTable$family, data.project.samples)
genus <- ChooseProjectSamples(totalTable$genus, data.project.samples)
species <- ChooseProjectSamples(totalTable$species, data.project.samples)

# todo: after applying:
#
# > length(colnames(species.ee))
# [1] 1824
# > length(unique(colnames(species.ee)))
# [1] 1703
colnames(species) <- make.short.name(colnames(species)) 
##

family.ee <- ChooseProjectSamples(family, meta.list$data.samples.paired.ee)
genus.ee <- ChooseProjectSamples(genus, meta.list$data.samples.paired.ee)
species.ee <- ChooseProjectSamples(species, meta.list$data.samples.paired.ee)

family.ka <- ChooseProjectSamples(family, meta.list$data.samples.paired.ka)
genus.ka <- ChooseProjectSamples(genus, meta.list$data.samples.paired.ka)
species.ka <- ChooseProjectSamples(species, meta.list$data.samples.paired.ka)

ee.features <- ChooseFeaturesEE(family.ee, genus.ee, species.ee)
ka.features <- ChooseFeaturesKA(family.ka, genus.ka, species.ka)

# bacteria top for EE group
outdir <- "/home/anna/metagenome/16s_antib/out_ee_ka/"

top.info.before <- GetTopInfo(species[rownames(species) %in% meta.list$data.samples.paired.ee$before,], "before")
top.info.after <- GetTopInfo(species[rownames(species) %in% meta.list$data.samples.paired.ee$after,], "after")

write.table(top.info.before$pool.top, 
            file = paste0(outdir, "ee_top_before.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

write.table(top.info.after$pool.top, 
            file = paste0(outdir, "ee_top_after.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

# bacteria top for KA group
top.info.before <- GetTopInfo(species[rownames(species) %in% meta.list$data.samples.paired.ka$before,], "before")
top.info.after <- GetTopInfo(species[rownames(species) %in% meta.list$data.samples.paired.ka$after,], "after")

write.table(top.info.before$pool.top, 
            file = paste0(outdir, "ka_top_before.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

write.table(top.info.after$pool.top, 
            file = paste0(outdir, "ka_top_after.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

#Top species in EE KA groups
# todo: how does it work? what does 85 mean?
top.species.ee.1 <- chooseTOPfeature(ee.features$species.ee.1, 85)
top.species.ee.2 <- chooseTOPfeature(ee.features$species.ee.2, 85)

top.species.ka.1 <- chooseTOPfeature(ka.features$species.ka.1, 85)
top.species.ka.2 <- chooseTOPfeature(ka.features$species.ka.2, 85)

top.species.ee <-list(top.species.ee.1=top.species.ee.1, top.species.ee.2=top.species.ee.2)
top.species.ka <-list(top.species.ee.1=top.species.ka.1, top.species.ee.2=top.species.ka.2)


mainDir <- "/home/anna/metagenome/16s_antib"
subDir <- "out_ee_ka"

if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
  
}

mainDir <- "/home/anna/metagenome/16s_antib/out_ee_ka"
subDir <- "heatmap"

if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
  
}


#Make HeatMap
HeatMap(top.species.ee.1, 'out_ee_ka/heatmap/', 'species_EE_1')
HeatMap(top.species.ee.2, 'out_ee_ka/heatmap/', 'species_EE_2')
HeatMap(top.species.ka.1, 'out_ee_ka/heatmap/', 'species_KA_1')
HeatMap(top.species.ka.2, 'out_ee_ka/heatmap/', 'species_KA_2')

#Print Mean and Stdev for group of samples 

PrintMeanAndStd(top.species.ee, min.trsh = 0.8, "EE" )
PrintMeanAndStd(top.species.ee, min.trsh = 0.8, "KA" )


# pool
pool_directory <- paste0("out_ee_ka/", "pool/")
if (!dir.exists(pool_directory)) {
  dir.create(pool_directory, showWarnings = F)
}

GroupPool(species,
          meta.data.g1 = meta.list$data.groups.ee[group_id == "EE_1"]$sample_id,
          meta.data.g2 = meta.list$data.groups.ee[group_id == "EE_2"]$sample_id,
          group.name = "EE")

GroupPool(species,
          meta.data.g1 = meta.list$data.groups.ka[group_id == "KA_1"]$sample_id,
          meta.data.g2 = meta.list$data.groups.ka[group_id == "KA_2"]$sample_id,
          group.name = "KA")

# PI-CRUST: prepare data
pred_table <- LoadPicrust(projects.directory.path='/home/anna/metagenome/16s_study/data/qiime/Case', 
                          project.ids=meta.list$data.groups.ee.ka$sample_id)
sample_cols <- 2:ncol(pred_table)
pred_table <- pred_table[, (sample_cols):=lapply(.SD, function(x) {x*100.0/sum(x)}), .SDcols=sample_cols]

#vitamin metabolism
outdir <- '/home/anna/metagenome/16s_study/out_ee_ka/PI-CRUST'
vitamin_table <- fread(pathway$kegg$kegg.vitamins, col.names=c('name_path', 'kegg_ko', 'vitamin'))
kegg_table <- fread(pathway$kegg$kegg.table, col.names=c('kegg_ko', 'kegg_id'))
kegg_id_vitamins <- merge(kegg_table, vitamin_table[, c('kegg_ko', 'vitamin'), with=F], by='kegg_ko')

PiCrustVitaminGroup(meta.1 = meta.list$data.groups.ee[group_id=="EE_1"]$sample_id,
             meta.2 = meta.list$data.groups.ee[group_id=="EE_2"]$sample_id,
             plot.name.1 = "EE_1_before",
             plot.name.2 = "EE_1_after",
             group.name = "EE")

PiCrustVitaminGroup(meta.1 = meta.list$data.groups.ka[group_id=="KA_1"]$sample_id,
             meta.2 = meta.list$data.groups.ka[group_id=="KA_2"]$sample_id,
             plot.name.1 = "KA_1_before",
             plot.name.2 = "KA_1_after",
             group.name = "KA")

#acids metabolism#
butyrate_table <- suppressWarnings(fread(pathway$kegg$kegg.acids, sep=',', 
                                         header=T, col.names=c('intermediate', 'kegg_id', 'name')))

pred_table <- LoadPicrust(projects.directory.path='/home/anna/metagenome/16s_study/data/qiime/Case', 
                          project.ids=meta.list$data.groups.ee.ka$sample_id)
sample_cols <- 2:ncol(pred_table)
pred_table <- pred_table[, (sample_cols):=lapply(.SD, function(x) {x*100.0/sum(x)}), .SDcols=sample_cols]

PiCrustButyrateGroup(meta.data.1 = meta.list$data.samples.paired.ee$before,
                     meta.data.2 = meta.list$data.samples.paired.ee$after,
                     group.name = "EE",
                     pred_data=pred_table, butyrate=butyrate_table,
                     plot.name.1 = "EE_1_before",
                     plot.name.2 = "EE_1_after")

PiCrustButyrateGroup(meta.data.1 = meta.list$data.samples.paired.ka$before,
                     meta.data.2 = meta.list$data.samples.paired.ka$after,
                     group.name = "KA",
                     pred_data=pred_table, butyrate=butyrate_table,
                     plot.name.1 = "KA_1_before",
                     plot.name.2 = "KA_1_after")

# alpha diversity
totalDiversity <- LoadAlphaDiv('/home/anna/metagenome/16s_study/data_ee_ka/qiime/Case')
alpha.div.ee1 <- ChooseAlphaSamples(totalDiversity, meta.list$data.groups.ee[group_id=='EE_1'])
alpha.div.ee2 <- ChooseAlphaSamples(totalDiversity, meta.list$data.groups.ee[group_id=='EE_2'])

alpha.div.ee1$group <- 1
alpha.div.ee2$group <- 2
alpha.div.ee <- rbind(alpha.div.ee1, alpha.div.ee2)

alpha.div.ka1 <- ChooseAlphaSamples(totalDiversity, meta.list$data.groups.ka[group_id=='KA_1'])
alpha.div.ka2 <- ChooseAlphaSamples(totalDiversity, meta.list$data.groups.ka[group_id=='KA_2'])
alpha.div.ka1$group <- 1
alpha.div.ka2$group <- 2
alpha.div.ka <- rbind(alpha.div.ka1, alpha.div.ka2)

MakeAlphaBoxPlot(alpha.div.ee, 'EE_alpha_diversity')
MakeAlphaBoxPlot(alpha.div.ka, 'KA_alpha_diversity')

#Statistics
feature.TOPlist.ee <- lapply(ee.features, chooseTOPfeature, perc=85)
feature.TOPlist.ka <- lapply(ka.features, chooseTOPfeature, perc=85)

family.ee.TOP <- chooseTOPfeature(family.ee, perc = 85)
family.ka.TOP <- chooseTOPfeature(family.ka, perc = 85)

genus.ee.TOP <- chooseTOPfeature(genus.ee, perc = 85)
genus.ka.TOP <- chooseTOPfeature(genus.ee, perc = 85)

species.ee.TOP <- chooseTOPfeature(species.ee, perc = 85)
species.ka.TOP <- chooseTOPfeature(species.ka, perc = 85)

outdir_stat = '/home/anna/metagenome/16s_study/out_ee_ka/stat/'

#KA family
StatResult(family.ka.TOP, meta.list$data.groups.ka[group_id=="KA_1"]$sample_id,
           meta.list$data.groups.ka[group_id=="KA_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "KA_family_wilcox", test = FALSE)

StatResult(family.ka.TOP, meta.list$data.groups.ka[group_id=="KA_1"]$sample_id,
           meta.list$data.groups.ka[group_id=="KA_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "KA_family_wilcox", test = TRUE)
#KA genus
StatResult(genus.ka.TOP, meta.list$data.groups.ka[group_id=="KA_1"]$sample_id,
           meta.list$data.groups.ka[group_id=="KA_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "KA_genus_ttest", test = FALSE)

StatResult(genus.ka.TOP, meta.list$data.groups.ka[group_id=="KA_1"]$sample_id,
           meta.list$data.groups.ka[group_id=="KA_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "KA_genus_wilcox", test = TRUE)
#KA species
StatResult(species.ka.TOP, meta.list$data.groups.ka[group_id=="KA_1"]$sample_id,
           meta.list$data.groups.ka[group_id=="KA_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "KA_species_ttest", test = FALSE)

StatResult(species.ka.TOP, meta.list$data.groups.ka[group_id=="KA_1"]$sample_id,
           meta.list$data.groups.ka[group_id=="KA_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "KA_species_wilcox", test = TRUE)


#EE family
StatResult(family.ee.TOP, meta.list$data.groups.ee[group_id=="EE_1"]$sample_id,
           meta.list$data.groups.ee[group_id=="EE_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "EE_ttest", test = FALSE)

StatResult(family.ee.TOP, meta.list$data.groups.ee[group_id=="EE_1"]$sample_id,
           meta.list$data.groups.ee[group_id=="EE_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "EE_wilcox", test = TRUE)
#EE genus
StatResult(genus.ee.TOP, meta.list$data.groups.ee[group_id=="EE_1"]$sample_id,
           meta.list$data.groups.ee[group_id=="EE_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "EE_ttest", test = FALSE)

StatResult(genus.ee.TOP, meta.list$data.groups.ee[group_id=="EE_1"]$sample_id,
           meta.list$data.groups.ee[group_id=="EE_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "EE_wilcox", test = TRUE)
#EE genus
StatResult(species.ee.TOP, meta.list$data.groups.ee[group_id=="EE_1"]$sample_id,
           meta.list$data.groups.ee[group_id=="EE_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "EE_ttest", test = FALSE)

StatResult(species.ee.TOP, meta.list$data.groups.ee[group_id=="EE_1"]$sample_id,
           meta.list$data.groups.ee[group_id=="EE_2"]$sample_id, maxpv = 0.05,
           pairedt = TRUE, outdir_stat, nameplot = "EE_wilcox", test = TRUE)



##################
# enterotype clusters analysis
message(paste(Sys.time(), "cluster analysis start"))
out_directory1 = '/home/anna/metagenome/16s_study/out_ee_ka/'

png(paste0(out_directory1, 'EE_KA.png'), width=800, height=800)
clusters.and.plot <- GetClustersAndPlot(abundance.matrix = species, data.samples = meta.list$data.samples.ee.ka.pair)
PlotBeforeAfterClusters(clusters.and.plot$bca.result, clusters.and.plot$data.sample.cluster, '', '')
invisible(dev.off())

message(paste(Sys.time(), "cluster analysis finish"))
