rm(list=ls(all=TRUE)); gc()

current_directory <- "/home/anna/metagenome/16s_study"
#current_directory <- "~/kids/16s_study"
#current_directory <- "/projects/16s_study"
#current_directory <- "/projects/16s_study"
#current_directory <- "~/do/16s_study"

setwd(current_directory)

list.of.packages <- c("stringr", 'ecodist', 'data.table', 'ape',
                      'ggplot2', 'scales', 'MASS', 'stringr', 
                      'matrixStats', 'gridExtra', 'grid', 
                      'reshape', 'gplots', 'MASS', 'futile.logger',
                      'lint', 'asbio', 'igraph', 'vegan', 'clusterSim')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")

# IT'S FORBIDDEN TO ADD FUNCTION TO FUNCTIONS.R!
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
source("lib/disease_facts.R")

library(data.table )
library(reshape)
library(plyr)
library(RColorBrewer)

message(paste(Sys.time(), "analysis start"))
path.to.config <- file.path(current_directory, 'src', 'pathway_ohmygut.ini')
pathway <- ReadIni(path.to.config) 

#loading meta
data.samples <- fread(pathway$meta.tables$meta.pairedsamples, colClasses=rep("character", 2))
data.groups <- fread(pathway$meta.tables$meta.samplesgroups, colClasses=rep("character", 3))
data.samples <- GroupSamples(data.samples, data.groups)
colnames(data.samples) <- c('before', 'after')

table.projects.by.sample <- fread(pathway$meta.tables$meta.readsets, sep=',', colClasses='character', 
                                  select=c(1,2), col.names = c('sample', 'run_id'))
table.projects.folders <- fread(pathway$meta.tables$meta.runs, sep='\t', colClasses='character', 
                                select = c(1, 4), col.names = c('run_id', 'folder'))
table.projects.by.sample <- merge(table.projects.by.sample, table.projects.folders, by='run_id')

totalTable <- LoadProjectsBySample(projects.directory.path=pathway$main$project.folder,
                                   table.projects = table.projects.by.sample,
                                   file.name.family = pathway$main$name.family,
                                   file.name.genus = pathway$main$name.genus,
                                   file.name.species = pathway$main$name.species,
                                   folder = pathway$main$inner.project.folder)



family <- totalTable$family
genus <- totalTable$genus
species <- totalTable$species


colnames(family) <- make.short.name(colnames(family))
family <- aggregate.by.name(family)
colnames(genus) <- make.short.name(colnames(genus))
genus <- aggregate.by.name(genus)
colnames(species) <- make.short.name(colnames(species))
species <- aggregate.by.name(species)

#only paired succesful reads from project
data.samples <- data.samples[(data.samples$before %in% rownames(species)) & (data.samples$after %in% rownames(species))]
data.groups <- data.groups[(data.groups$sample_id %in% c(data.samples$before, data.samples$after))]
species <- species[rownames(species) %in% data.groups$sample_id,]
genus <- genus[rownames(genus) %in% data.groups$sample_id,]
family <- family[rownames(family) %in% data.groups$sample_id,]

###
# create sample directories
out.directory <- pathway$outputdir$outdir
by.sample.directory <- file.path(out.directory, "by_sample")
if (!dir.exists(out.directory)) {
  dir.create(out.directory, showWarnings = F)
}
if (!dir.exists(by.sample.directory)) {
  dir.create(by.sample.directory, showWarnings = F)
}

data.samples$dir <- file.path(by.sample.directory, paste0(data.samples$before, "__", data.samples$after))
sapply(data.samples$dir, function(x) if(!dir.exists(x)){ dir.create(x, showWarnings = F) })

MIN_CHANGE_RATIO <- 1.5 # minimum fold change (e.g. increase OR decrease )
MIN_MAX_PERC <- 0.01 # minimum threshold for maximum of taxon rel. abundance across 2 time points

########################
########################
#START ANALYSIS
########################
########################

# pool
pool_directory <- file.path(pathway$outputdir$outdir, "pool/")
if (!dir.exists(pool_directory)) {
  dir.create(pool_directory, showWarnings = F)
}
pool.before <- GetPool(species[rownames(species) %in% data.samples$before,], MIN_MAX_PERC)
pool.after <- GetPool(species[rownames(species) %in% data.samples$after,], MIN_MAX_PERC)

write.table(pool.before, 
            file = paste0(pool_directory, "all-pool-before.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

write.table(pool.after, 
            file = paste0(pool_directory, "all-pool-after.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

##################
# bacteria top & top pool
top.info.before <- GetTopInfo(species[rownames(species) %in% data.samples$before,], "before")
top.info.after <- GetTopInfo(species[rownames(species) %in% data.samples$after,], "after")
top.table <- cbind(top.info.before$top.by.sample.table, top.info.before$top.by.sample.table)

write.table(top.info.before$pool.top, 
            file = file.path(pool_directory, "pool-top-before.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

write.table(top.info.after$pool.top, 
            file = file.path(pool_directory, "pool-top-after.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

write.table(top.table, 
            file = file.path(pool_directory, "toptable-species.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

##################
# most increased/decreased bacterias
tmp <- data.samples[,{
  bugs <- get.changed(species, before, after, MIN_MAX_PERC, MIN_CHANGE_RATIO)
  write.table(bugs$taxa_decreased, file.path(dir, paste0("top_decreased_", before, "_", after, ".csv")),
              sep = "\t", row.names = F, quote = F, col.names = T)
  write.table(bugs$taxa_increased, file.path(dir, paste0("top_increased_", before, "_", after, ".csv")),
              sep = "\t", row.names = F, quote = F, col.names = T)
}, by = dir]

##################
# top histograms
data.samples[,{
  top.plot.before <- GetTopPlot(top.info.before$top.by.sample[[.SD$before]])
  ggsave(plot = top.plot.before, 
         filename = paste0(dir, "/species_top_before.png"), 
         width=4, height=2)
  top.after <- GetTopPlot(top.info.after$top.by.sample[[.SD$after]])
  ggsave(plot = top.after, 
         filename = paste0(dir, "/species_top_after.png"), 
         width=4, height=2)
  message(paste("writing top plot to:", dir))
}, by = dir]

##################
# rank abundance curve
top.by.sample.dt.before <- GetRankDataTable(top.info.before$top.by.sample)
base.plot.before <- GetRankAbundenceCurveBasePlot(top.by.sample.dt.before)
# example
GetHighlitedRankPlot(base.plot.before, 
                     top.by.sample.dt.before,
                     sample.id = "1")

top.by.sample.dt.after <- GetRankDataTable(top.info.after$top.by.sample)
base.plot.after <- GetRankAbundenceCurveBasePlot(top.by.sample.dt.after)

data.samples[,{
  rac_before <- GetHighlitedRankPlot(base.plot.before, 
                                     top.by.sample.dt.before,
                                     sample.id = .SD$before)
  ggsave(plot = rac_before, 
         filename = paste0(dir, "/species_rac_before.png"),
         width=3, height=2)
  rac_after <- GetHighlitedRankPlot(base.plot.after, 
                                    top.by.sample.dt.after,
                                    sample.id = .SD$after)
  ggsave(plot = rac_after, 
         filename = paste0(dir, "/species_rac_after.png"),
         width=3, height=2)
  
  message(paste("writing rank ab. curve to:", dir))
}, by = dir]

##################
# enterotype clusters analysis
message(paste(Sys.time(), "cluster analysis start"))

clusters <- GetClusters(abundance.matrix = species, data.samples = data.samples)

# first just plot and save clusters
png(file.path(out.directory, 'clusters.png'), width=1200, height=800)
PlotClusters(clusters$bca.result, clusters$data.cluster, clusters$cluster.perc$driver)
invisible(dev.off())

# then plot clusters for each sample
data.samples[,{
  png(file.path(dir, 'clusters.png'), width=1200, height=800)
  PlotClusters(clusters$bca.result, clusters$data.cluster, clusters$cluster.perc$driver)
  PlotBeforeAfterClusters(clusters$bca.result, clusters$data.sample.cluster, .SD$before, .SD$after)
  invisible(dev.off())
  message(paste('clusters plot for', .SD$before, .SD$after, 'is ready'))
}, by = dir]

# writing clusters info
write.table(clusters$cluster.perc, file.path(out.directory, "cluster_perc.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)
write.table(clusters$cluster.moves, file.path(out.directory, "cluster_moves.csv"), 
            sep = "\t", row.names = F, quote = F, col.names = T)

message(paste(Sys.time(), "cluster analysis finish"))

##################
# CO-OCCURANCE: prepare data
message(paste(Sys.time(), "coocurance analysis start"))
layout.path=pathway$cooc$cooc.coords
INTERACTIVE_MODE_SAVE_LAYOUT <- F
NORM_WEIGHTS <- T
features.genus <- copy(genus)

# CO-OCCURENCE: calculate 
cooccurance.all <- get.coocurance(features.genus, min.genera.percent=MIN_MAX_PERC, min.abs.corr=0.4, max.abs.corr=1)
# CO-OCCURENCE: plot total
png(file.path(out.directory, 'total_cooccur.png'), width=2500, height=1800)
plot.graph(cooccurance.all, layout.path = layout.path, 
           save.layout=INTERACTIVE_MODE_SAVE_LAYOUT, alpha.v=0.8, alpha.e=0.3, 
           scale.v = function(x) 15*x^0.3, scale.e = function(x) 80*(abs(x)^4))
###################STOP IF INTERACTIVE_MODE_SAVE_LAYOUT
invisible(dev.off())

if (NORM_WEIGHTS) {
  coef <- 0.02/cooccurance.all$vertex$weight
  names(coef) <- cooccurance.all$vertex$name
  cooccurance.all$vertex$weight <- 0.02 * cooccurance.all$vertex$weight/cooccurance.all$vertex$weight
} else {
  coef <- rep(1, nrow(cooccurance.all$vertex))
  names(coef) <- cooccurance.all$vertex$name
}

# CO-OCCURANCE: plot for sample
# takes a few minutes
data.samples[, {
  png(file.path(dir, paste0(before, '_cooccur' , NORM_WEIGHTS, '.png')),  width=2500, height=1800)
  plot.for.each.sample(cooccurance.all, layout.path, features.genus[before, names(coef)] * coef, alpha.v=0.5, alpha.e=0.2, 
                       scale.v = function(x) 15*x^0.3, scale.e = function(x) 80*(abs(x)^4))
  invisible(dev.off())
  png(file.path(dir, paste0(after, '_cooccur' , NORM_WEIGHTS, '.png')),  width=2500, height=1800)
  plot.for.each.sample(cooccurance.all, layout.path, features.genus[after, names(coef)] * coef, alpha.v=0.5, alpha.e=0.2, 
                       scale.v = function(x) 15*x^0.3, scale.e = function(x) 80*(abs(x)^4))
  invisible(dev.off())
  cat("Create graph for ", dir, "\n")
}, by = 1:nrow(data.samples)]
message(paste(Sys.time(), "coocurance analysis finish"))

##################
# PI-CRUST: prepare data
pred_table <- LoadPicrustBySample(projects.directory.path=pathway$main$project.folder,
                                  project.ids=data.groups$sample_id,
                                  table.projects = table.projects.by.sample)
sample_cols <- 2:ncol(pred_table)
pred_table <- pred_table[, (sample_cols):=lapply(.SD, function(x) {x*100.0/sum(x)}), .SDcols=sample_cols]

#vitamin metabolism
vitamin_table <- fread(pathway$kegg$kegg.vitamins, col.names=c('name_path', 'kegg_ko', 'vitamin'))
kegg_table <- fread(pathway$kegg$kegg.table, col.names=c('kegg_ko', 'kegg_id'))
picrust_analysis_vitamins(pred_data=pred_table, 
                          quiime.data=data.samples, 
                          vitamins=vitamin_table, 
                          kegg_ko=kegg_table)

#acids metabolism
butyrate_table <- suppressWarnings(fread(pathway$kegg$kegg.acids, sep=',', header=T, col.names=c('intermediate', 'kegg_id', 'name')))
picrust_analysis_butyrate(pred_data=pred_table, 
                          quiime.data=data.samples, 
                          butyrate=butyrate_table)

message(paste(Sys.time(), "analysis finish"))

##################
# Alpha diversity

totalDiversity <- LoadAlphaDiv(file.path(current_directory, 'data_2016_04_28/qiime/2015_04_28_gg_base'))

alpha.div.project <- ChooseAlphaSamples(totalDiversity, data.groups)

alpha.div.project.g1 <- alpha.div.project[alpha.div.project$sample %in% data.samples$before]
alpha.div.project.g2 <- alpha.div.project[alpha.div.project$sample %in% data.samples$after]
alpha.div.project.g1$group <- 1
alpha.div.project.g2$group <- 2
alpha.div.project <- rbind(alpha.div.project.g1, alpha.div.project.g2)
alpha.div.project$diversity_scaled <- (alpha.div.project$diversity / max(alpha.div.project$diversity)) * 10
# group is a factor, not a continuos variable!!
alpha.div.project$group <- factor(alpha.div.project$group, labels=c("before", "after"))

base.plot <- GetAlphaDivBasePlot(alpha.div.project)

data.samples[, {
  p_id <- GetAlphaDivSamplePlot(base.plot, data=alpha.div.project, sample.name=before )
  p_id <- GetAlphaDivSamplePlot(p_id, data=alpha.div.project, sample.name=after)
  ggsave(file.path(dir, paste0(before, '_', after, '_alpha.div.png')), device='png', width=10, height=5, plot = p_id)
  cat("Create graph for ", dir, "\n")
}, by = 1:nrow(data.samples)]

##################
# NeededSampleSize (in progress)

genus.bef <- genus[rownames(genus) %in% data.samples$after,]
genus.aft <- genus[rownames(genus) %in% data.samples$before,]

#собираем малопредставленные виды в один столбец, умножаем сначала данные на константу, 
#чтобы все ненулевые элементы матрицы были больше 1
#ToDo: сделать сранение всех групп из набора между собой
filt.before <- Data.filter(genus.bef*6000, "data", 1, 30)
filt.after <- Data.filter(genus.aft*6000, "data", 1, 30)

data <- merge.features(filt.before,filt.after, 0)
group.before <- data$dataCase
group.after <- data$dataCntrl

power.curve(group.before, group.after, 1 ,2000)
flog.info("END")

##################
# GET DISEASE FACTS:
fact.table = fread('text_mining_data/bact_disease_facts.csv', sep='\t')
all.bact <- fread('text_mining_data/all_bact_catalog.csv', sep=',')
tmp <- write.facts(tax.matrix=species,
                   dir.sample.id=data.samples,
                   fact.table=fact.table,
                   bact.catalog=all.bact)
tmp <- write.facts(tax.matrix=genus,
                   dir.sample.id=data.samples,
                   fact.table=fact.table,
                   bact.catalog=all.bact,
                   cut.names=T,
                   append = T)
tmp <- write.facts(tax.matrix=family,
                   dir.sample.id=data.samples,
                   fact.table=fact.table,
                   bact.catalog=all.bact,
                   cut.names=T,
                   append = T)








