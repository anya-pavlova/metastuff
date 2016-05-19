library(vegan)

top.species.ee <- chooseTOPfeature(species.ee, 85) # выбираются не нулевые бактерии для образцов ант.1.до+ант.1.после
top.species.ka <- chooseTOPfeature(species.ka, 85) # выбираются не нулевые бактерии для образцов ант.2.до+ант.2.после

s <- UniteMatrices(top.species.ee, top.species.ka) # объединяю их в одну матрицу
s <- t(s) # транспонирую, так как adonis принимает факты, по которым смотреть различия, по строкам (как это правильно говорить? Я очень криво выражаюсь, интересно бы знать нормальные названия и может более правильную логику) 
fac <- as.character(rownames(s)) # формирую фактор из названий бактерий
per <- adonis(s~fac, permutations = 9999, method = "bray") #запускаю сам тест


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
