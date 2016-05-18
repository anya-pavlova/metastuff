LoadProjects <- function(projects.directory.path, file.name.family, file.name.genus, file.name.species, folder = ''){
  if (!(folder == '')) file.name.family <- file.path(folder, file.name.family)
  if (!(folder == '')) file.name.genus <- file.path(folder, file.name.genus)
  if (!(folder == '')) file.name.species <- file.path(folder, file.name.species)
  dirs <- list.dirs(projects.directory.path, recursive=FALSE, full.names = TRUE)
  matrix.family <- NA 
  matrix.genus <- NA
  matrix.species <- NA
  for (dir in dirs){
    family <- file.path(dir, file.name.family)
    genus <- file.path(dir, file.name.genus)
    species <- file.path(dir, file.name.species)
    message(family)
    message(genus)
    message(species)
    matrix.family <- UniteMatrices(ReadQiimeSumFeats(family), matrix.family)
    matrix.genus <- UniteMatrices(ReadQiimeSumFeats(genus), matrix.genus)
    matrix.species <- UniteMatrices(ReadQiimeSumFeats(species), matrix.species)
  }
  return(list(family=matrix.family, genus=matrix.genus, species=matrix.species))
}

LoadProjectsBySample <- function(projects.directory.path, table.projects, 
                                 file.name.family, file.name.genus, file.name.species, folder = ''){
  if (!(folder == '')) file.name.family <- file.path(folder, file.name.family)
  if (!(folder == '')) file.name.genus <- file.path(folder, file.name.genus)
  if (!(folder == '')) file.name.species <- file.path(folder, file.name.species)
  families <- get.quiime.for.taxon(projects.directory.path, table.projects, file.name.family)
  genus <- get.quiime.for.taxon(projects.directory.path, table.projects, file.name.genus)
  species <- get.quiime.for.taxon(projects.directory.path, table.projects, file.name.species)
  return(list('family'=families, 'genus'=genus, 'species'=species))
}

get.quiime.for.taxon <- function(projects.directory.path, table.projects, file.name) {
  taxons <- alply(unique(table.projects$folder), 1, function (x) {
    taxon <- file.path(projects.directory.path, x, file.name)
    samples <- table.projects[folder==x]$sample
    message(taxon)
    data.taxon <- fread(taxon, header=T, sep = '\t', select=c('Taxon', samples))
  })
  #merge
  merged.taxons = taxons[[1]]
  for (taxon in taxons[-1]) {
    merged.taxons <- merge(merged.taxons, taxon, by = 'Taxon', all=T)
  }
  #make matrix
  matrix.taxons <- t(as.matrix(merged.taxons[, !'Taxon', with=F]))
  colnames(matrix.taxons) <- merged.taxons$Taxon
  matrix.taxons[is.na(matrix.taxons)]=0.0
  return(matrix.taxons)
}


#meta table for EE KA project
MetaLoading <- function(pathway.pair, patway.data){
  data.samples.all <- fread(pathway.pair, colClasses=rep("character", 2))
  data.groups.all <- fread(patway.data, colClasses=rep("character", 3))
  colnames(data.samples.all) <- c('before', 'after')
  
  data.groups.ee <- rbind(data.groups.all[group_id == "EE_1"], data.groups.all[group_id == "EE_2"])
  data.groups.ka <- rbind(data.groups.all[group_id == "KA_1"], data.groups.all[group_id == "KA_2"])
  data.groups.ee.ka <- rbind(data.groups.ee, data.groups.ka)
  data.samples.ee.ka.pair <- data.samples.all[data.samples.all$before %in% data.groups.ee.ka$sample_id]
  
  just.ee <- data.groups.all[group_id == "EE_1"]
  just.ka <- data.groups.all[group_id == "KA_1"]
  
  data.samples.paired.ee <- data.samples.ee.ka.pair[data.samples.ee.ka.pair$before %in% just.ee$sample_id]
  data.samples.paired.ka <- data.samples.ee.ka.pair[data.samples.ee.ka.pair$before %in% just.ka$sample_id]
  
  list(data.groups.ee = data.groups.ee, data.groups.ka = data.groups.ka, data.groups.ee.ka=data.groups.ee.ka,
       data.samples.ee.ka.pair = data.samples.ee.ka.pair,
       data.samples.paired.ee = data.samples.paired.ee, data.samples.paired.ka = data.samples.paired.ka,
       data.samples.ee.ka.pair = data.samples.ee.ka.pair
  )
}