library(vegan)

top.species.ee <- chooseTOPfeature(species.ee, 85)
top.species.ka <- chooseTOPfeature(species.ka, 85)

a1 <- head(top.species.ee,1)
a2 <- head(top.species.ka,3)

t.aa <- UniteMatrices(a1,a2)
#t.aa <- UniteMatrices(top.species.ee, top.species.ka)

t.aa <- t.aa[,c(1:3)]

t.aa.2 <- t(t.aa)

fac <- as.character(rownames(t.aa.2))
per <- adonis(t.aa.2~fac, permutations = 9999, method = "bray")



vegdist(t.aa, method = 'bray')


library(vegan)
library(MASS)
data(varespec)
vare.dis <- vegdist(varespec)
hh <- vegan::vegdist(t.aa)
t.aa[1:2,1:2]
rowSums(t.aa)

bcdist(t.aa)
dist.all<-bcdist(bact.for.samp.id.all)
dist.all.m <- as.matrix(dist.all)



#test
library(vegan)

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


betad <- betadiver(dune, "z")
adonis(betad ~ Management+A1, dune.env, perm=200)

ad <- adonis(betad ~ Management, dune.env, perm=200)
ad
str(ad)
ad$aov.tab$`Pr(>F)`[1]
