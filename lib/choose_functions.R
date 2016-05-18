#Choose functions
ChooseProjectSamples <- function(feature, data.samples)
{
  all.samples.id <- c(data.samples$before, data.samples$after)
  bact.for.all.fam <- feature[all.samples.id, ]
  bact.for.all.fam
}

ChooseGroupBefore <- function(feature, data.samples)
{
  all.samples.id <- data.samples$before
  bact.for.all.fam <- feature[all.samples.id, ]
  bact.for.all.fam
}

ChooseGroupAfter <- function(feature, data.samples)
{
  all.samples.id <- data.samples$after
  bact.for.all.fam <- feature[all.samples.id, ]
  bact.for.all.fam
}

ChooseFeaturesEE <- function(family.ee, genus.ee, species.ee)
{
  family.ee.1 <- ChooseGroupBefore(family.ee, meta.list$data.samples.paired.ee)
  genus.ee.1 <- ChooseGroupBefore(genus.ee, meta.list$data.samples.paired.ee)
  species.ee.1 <- ChooseGroupBefore(species.ee, meta.list$data.samples.paired.ee)
  
  family.ee.2 <- ChooseGroupAfter(family.ee, meta.list$data.samples.paired.ee)
  genus.ee.2 <- ChooseGroupAfter(genus.ee, meta.list$data.samples.paired.ee)
  species.ee.2 <- ChooseGroupAfter(species.ee, meta.list$data.samples.paired.ee)
  list(family.ee.1=family.ee.1, genus.ee.1=genus.ee.1, species.ee.1=species.ee.1,
       family.ee.2=family.ee.2, genus.ee.2=genus.ee.2, species.ee.2=species.ee.2)
}

ChooseFeaturesKA <- function(family.ka, genus.ka, species.ka)
{
  family.ka.1 <- ChooseGroupBefore(family.ka, meta.list$data.samples.paired.ka)
  genus.ka.1 <- ChooseGroupBefore(genus.ka, meta.list$data.samples.paired.ka)
  species.ka.1 <- ChooseGroupBefore(species.ka, meta.list$data.samples.paired.ka)
  
  family.ka.2 <- ChooseGroupAfter(family.ka, meta.list$data.samples.paired.ka)
  genus.ka.2 <- ChooseGroupAfter(genus.ka, meta.list$data.samples.paired.ka)
  species.ka.2 <- ChooseGroupAfter(species.ka, meta.list$data.samples.paired.ka)
  list(family.ka.1=family.ka.1, genus.ka.1=genus.ka.1, species.ka.1=species.ka.1,
       family.ka.2=family.ka.2, genus.ka.2=genus.ka.2, species.ka.2=species.ka.2)
}
