# random abundance matrix
v <- runif(n=80,min=0,max=1)
m <- matrix(v, nrow=4, ncol=20)
m <- m / apply(m,MARGIN=1,FUN=sum)
rownames(m) <- paste("sample_",c(1:4),sep="")
colnames(m) <- paste("bacteria_",c(1:20),sep="")

# abundance dataframe
df <- melt(m)
names(df) <- c("sample", "bacteria", "abundance")

#Plot1
ggplot(df,aes(x=sample,y=abundance, fill=bacteria))+
  geom_histogram(stat="identity",colour = "black")+
  scale_fill_manual(values=rainbow(n=20))+
  ylab("Relative abundance")+xlab("Samples")+
  theme_bw()

#Plot2
ggplot(df,aes(x=sample,y=abundance, fill=bacteria))+
  geom_histogram(stat="identity",colour = "black")+
  scale_fill_manual(values=rainbow(n=20))+
  ylab("Relative abundance")+xlab("Samples")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45))

