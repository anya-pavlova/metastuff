#
# bacterial co-occurence graph
#

# some parameters
MIN_GENERA_PERCENT <- 1 # drop low count genera
MIN_ALLOWED_ABS_CORR <- 0.4 # drop edges with lower absolute correlations
MAX_ALLOWED_ABS_CORR <- 1 #0.9

# create output directories
dir.create("out/cooccur")
dir.create("out/cooccur/per_sample/")
dir.create("out/cooccur/graphs/")
dir.create("out/cooccur/graphs/per_sample/")


# read in features

# this can be at any taxonomic level. shown for genus.
INPUT_FEATURES <- genus #projectTable$genus.g1 # from Anya's main.R

# shorten bac names
# TODO wrap it into a function
#colnames(g) <- paste("o_", unlist(data.frame(strsplit(colnames(g), "o_"))[2,]), sep="")
colnames(INPUT_FEATURES) <- paste("o_", unlist(data.frame(strsplit(colnames(INPUT_FEATURES), "o_"))[2,]), sep="")

g <- INPUT_FEATURES

# drop low count taxons
g <- g[,which(apply(g, 2, max) >= MIN_GENERA_PERCENT)]
meancov <- t(apply(g, 2, mean))

n_samples <- nrow(g)
n_feats <- ncol(g)

print(n_samples)
print(n_feats)

# can take only most prevalent microbes - when list is sorted by MAX abundance it is easy
#best_feats=c(1:24)
#g <- g[best_feats]
best_feats=c(1:ncol(g))

n_feats <- length(best_feats)

print("n_feats changed to")
print(n_feats)

# Spearman correlation 
corm <- cor(g, NULL, "everything", "spearman")
#heatmap.2(corm)

# create graph nodes
count <- 0
nx <- n_feats
nodes <- matrix(, nrow=nx*nx, ncol=5)
meancov_out <- matrix(, nrow=nx, ncol=2)
for(i in 1:(nx-1))
{
  for(j in (i+1):nx)
  {
    # drop extraordinary correlation values
    #if(abs(corm[i,j]) <= MAX_ALLOWED_ABS_CORR && abs(corm[i,j]) >= MIN_ALLOWED_ABS_CORR)
    
    # drop extraordinary and negative correlation values
    if(abs(corm[i,j]) <= MAX_ALLOWED_ABS_CORR && abs(corm[i,j]) >= MIN_ALLOWED_ABS_CORR && corm[i,j] >= 0)    
    {
      count<-count+1      
      nodes[count,] <- c(colnames(g)[i], colnames(g)[j], "cor_abs", abs(corm[i,j]), sign(corm[i,j]))
    }	
    else
    {	
      print("dropped")
      print(c(colnames(g)[i], colnames(g)[j], "cor", corm[i,j]))
    }
  }	
}

# prepare means
for(k in 1:nx)
{
  meancov_out[k,] <- c(colnames(g)[k], meancov[1, best_feats[k]])
}

# append singletons
a <- rbind(nodes[,1], nodes[,2])
b <- which(!(meancov_out[,1] %in% a))
for(j in b)
{
  count<-count+1      
  nodes[count,] <- c(meancov_out[j,1], meancov_out[j,1], "cor_abs", 1, 1)
}

# remove NAs
ind<-which(apply(nodes,1,function(x)all(is.na(x))))
if(length(ind) != 0)
{
  nodes<-nodes[-ind,]
}

# save the graphs as digits for external usage (e.g., in Cytoscape)
write(t(nodes), file = "out/cooccur/corr_nodes.txt", ncolumns=5, append = FALSE, sep = "\t")
write(t(meancov_out), file = "out/cooccur/meancov_out.txt", ncolumns=2, append = FALSE, sep = "\t")

# generate weighted graph for each sample
sample_weights <- g[,meancov_out[,1]]
#rownames(sample_weights) <- names(meancov_out)
#colnames(sample_weights) <- rownames(g)
#write(sample_weights, file = "out/sample_weights.txt", append = FALSE, sep = "\t")
write.table(round(sample_weights, 5), file = "out/cooccur/sample_weights.txt", sep="\t", row.names=TRUE, quote=FALSE, col.names = NA)
for(i in 1:n_samples)
{
  write.table(round(t(sample_weights), 5)[,i], file = paste("out/cooccur/per_sample/sw_", rownames(sample_weights)[i], ".txt", sep=""), sep="\t", row.names=TRUE, quote=FALSE, col.names = FALSE)
}


# get max
maxcov <- apply(INPUT_FEATURES, 2, max)[colnames(meancov)]


# use igraph for visualization
nodes2 <- data.frame(nodes[,c(1,2,4)], stringsAsFactors=F)
colnames(nodes2) <- c("bac1", "bac2", "cor_abs")

# remove isolated nodes for now
nodes2 <- nodes2[-which(nodes2[,1] == nodes2[,2]),]

nodes2$cor_abs <- as.numeric(nodes2$cor_abs)
net <- graph.data.frame(data.frame(nodes2), directed=F)
net <- set.vertex.attribute(net, "mean", value=meancov[,V(net)$name])
net <- set.vertex.attribute(net, "max", value=maxcov[V(net)$name])

#net <- simplify(net)
# randomize 
set.seed(123)
# experiments with layout
#net <- set.graph.attribute(net, "layout", layout.fruchterman.reingold.grid(net, repulserad=10000))
#net <- set.graph.attribute(net, "layout", layout.fruchterman.reingold.grid(net, niter=10000))
lay <- layout.fruchterman.reingold.grid(net, repulserad=10000)
lay <- layout.fruchterman.reingold.grid(net, niter=10000)
lay2 <- layout.norm(lay, -10, 10, -10, 10)
#lay <- layout.fruchterman.reingold.grid(net, niter=10000)
#lay <- layout.norm(lay, -0.1, 0.1, -0.1, 0.1)
net <- set.graph.attribute(net, "layout", lay)

par(mai=c(0,0,0,0)) 

# draw graph in a window
tkp <- tkplot(net, vertex.size=0.5*V(net)$max, edge.width=4*get.edge.attribute(net, name="cor_abs"), vertex.label.cex=0.9,
              edge.color="grey", vertex.color="red", vertex.label.color="black", vertex.label.family="Arial", vertex.frame.color="grey",
              canvas.width=1000, canvas.height=800)
tkplot.fit.to.screen(tkp)

########################################################
# User, interactively adjust the positions here!
# so that the graph looks nice, labels do not overlap, etc.
# Then the same vertex coordinates will be used to draw all graphs to files automatically.
########################################################

# get the adjuasted coords and draw them
coords <- tkplot.getcoords(tkp)
# SWITCH: 
# save...
save(coords, file="out/cooccur/manual_coords.RData")
# .. or, once precomputed nicely, just load them from cache
#load(file="out/cooccur/manual_coords.RData")

# save group graph
png(file=paste("out/cooccur/graphs/_gr_TOTAL.png", sep=""), width=1000, height=1000)
  plot.igraph(net, vertex.size=0.5*V(net)$max, edge.width=4*get.edge.attribute(net, name="cor_abs"), vertex.label.cex=0.9,
            edge.color="grey", vertex.color="red", vertex.label.color="black", vertex.frame.color="grey", 
            vertex.label.family="Arial", layout=coords, margin=0)
dev.off()

# save individual graphs
for(sid in rownames(INPUT_FEATURES))
{
  sampleNet <- net
  sampleNet <- set.vertex.attribute(sampleNet, "mean", value=unlist(INPUT_FEATURES[sid, V(net)$name]))
  V(sampleNet)$mean
  set.seed(108)

  png(file=paste("out/cooccur/graphs/per_sample/gr_", sid, ".png", sep=""), width=1000, height=1000)
  par(mai=c(0,0,0,0)) 
  plot.igraph(sampleNet, vertex.size=V(sampleNet)$mean, edge.width=3*get.edge.attribute(sampleNet, name="cor_abs"),vertex.label.cex=0.9,
              edge.color="grey", vertex.color="red", vertex.label.color="black", vertex.label.family="Arial", vertex.frame.color="grey",
              layout=coords, margin=0)    
  #plot.igraph(sampleNet, vertex.size=V(sampleNet)$mean) #, edge.width=3*get.edge.attribute(sampleNet, name="cor_abs"), vertex.label.cex=0.7)    
  dev.off()
}

