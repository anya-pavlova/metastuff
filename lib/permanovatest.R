library(vegan)

species.ee <- species[meta.list$data.groups.ee$sample_id,]
species.ka <- species[meta.list$data.groups.ka$sample_id,]

total.m <- UniteMatrices(species.ee, species.ka)
total.m.t <- t(total.m)
fac <- as.character(rownames(total.m.t))

per <- adonis(total.m.t~fac)
per <- adonis(rbind(dat1a,dat1b)~fac)

t.spe.ee <- t(as.data.table(species.ee, keep.rownames = T))
t.spe.ka <- t(as.data.table(species.ka, keep.rownames = T))

t.spe.ee[,1]
t.spe.ka[,1]

t.spe.ee <-  t.spe.ee[-1,]
t.spe.ka <-  t.spe.ka[-1,]

fac <- as.factor(c(rownames(t.spe.ee), rownames(t.spe.ka)))

kk <-UniteMatrices(species.ka, species.ee)



t.spe.ee <- t(species.ee)
t.spe.ka <- t(species.ka)



dat1a<-matrix(sample(c(0,1,1,1),200,replace=T),10,20)
dat1b<-matrix(sample(c(0,1,1,1),200,replace=T),10,20)

dat2<-matrix(sample(c(0,0,0,1),200,replace=T),10,20)

fac<-gl(2,10)
fac
dist11<-vegdist(rbind(dat1a,dat1b))
dist12<-vegdist(rbind(dat1a,dat2))
head(dist11)

an <- anova(betadisper(dist11,fac))
per <- adonis(rbind(dat1a,dat1b)~fac)

anova(betadisper(dist12,fac))
adonis(rbind(dat1a,dat2)~fac)

#test
a <- matrix(rep(4, 12), 2,6)
b <- matrix(rep(3, 12), 2,6)
fac1 <- (fac<-gl(2,2))
per <- adonis(rbind(dat1a,dat1b)~fac)



betad <- betadiver(dune, "z")
adonis(betad ~ Management+A1, dune.env, perm=200)

ad <- adonis(betad ~ Management, dune.env, perm=200)
ad
str(ad)
ad$aov.tab$`Pr(>F)`[1]