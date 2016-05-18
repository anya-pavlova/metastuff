library(igraph)
library(data.table)

get.coocurance <- function(input.features.matrix, min.genera.percent=1, min.abs.corr=0.4, max.abs.corr=1, verbose=T) {
  if (verbose) cat("Start calculating data for co-occarance graph...\n")
  
  features.matrix <- input.features.matrix
  features.matrix <- features.matrix[, apply(features.matrix, 2, max) >= min.genera.percent]
  n.feats <- ncol(features.matrix)
  
  cor.matrix <- cor(features.matrix, method = "spearman", use='everything')
  compare.cols <- function(i, j, namex, namey) {
    data.frame(name1=namex, name2=namey, value=cor.matrix[i, j], stringsAsFactors = F)
  }
  
  # MAKE CORRELATION 
  net <- ldply(1:(n.feats-1), function(i) {
    ldply((i+1):n.feats, function(j) compare.cols(i, j, colnames(features.matrix)[i], 
                                                  colnames(features.matrix)[j]))})
  net <- net[(abs(net$value)>=min.abs.corr) & (abs(net$value)<=max.abs.corr), ]
  
  # MAKE VERTEX ATTRIBUTES: weight, groups
  features.matrix <- features.matrix[, colnames(features.matrix) %in% c(net$name1, net$name2)]
  vertex <- ldply(1:ncol(features.matrix), function(x) {
    data.frame(name=colnames(features.matrix)[x], weight=mean(features.matrix[,x]))
  })
  if (verbose) cat("Finish calculating data for co-occarance graph...\n")
  return(list("vertex"= vertex, "net"=net))
}

add.alpha <- function(col, alpha=1){
  rgb.col = (col2rgb(col)/255)
  return(rgb(rgb.col['red',], rgb.col['green',], rgb.col['blue',], alpha=alpha))
}


##################################################
#visualization
##################################################
interactive.layout.plan <- function(net, path, scale.v = function(x) 10*x^0.5, koef=5) {
  lay <- layout_nicely(net, niter=5000)
  tkp = tkplot(net, layout = lay, vertex.size = koef*scale.v(V(net)$weight), canvas.width=1000, canvas.height=1000)
  readline(prompt="Press [enter] to continue")
  layout <- tkplot.getcoords(tkp)
  save(layout, file=path)
  tk_off()
}

add.alpha <- function(col, alpha=1){
  rgb.col = (col2rgb(col)/255)
  return(rgb(rgb.col['red',], rgb.col['green',], rgb.col['blue',], alpha=alpha))
}

fit.to.layout <- function(lay) {
  p.xlim = c(min(lay[,1]), max(lay[,1]))
  p.ylim = c(min(lay[,2]), max(lay[,2]))
  reshape_x=p.xlim[2]-p.xlim[1]
  reshape_y=p.ylim[2]-p.ylim[1]
  layout.out <- lay
  layout.out[,1] <- layout.out[,1]/reshape_x
  layout.out[,2] <- layout.out[,2]/reshape_y
  return(layout.out)
}

####################3
#vertex.vis may contain: name, weight, group, color, sec.group, frame.color
####################
plot.graph <- function(net.vis, layout.path=NULL, save.layout=F, alpha.v=0.7, alpha.e=0.2, 
                       scale.v = function(x) 10*x^0.5, scale.e = function(x) 20*(abs(x)^0.5), legend_x= 0.7, 
                       legend_y=0.1, cex_legend=2.4) {
  net <- graph.data.frame(net.vis$net, directed=F, vertices = net.vis$vertex)
  
  if (save.layout) {
    cat("Warning: use interactive mode for saving layout!!!")
    interactive.layout.plan(net, layout.path, scale.v = scale.v, koef = 5)
    return(0)
  } 
  
  # VERTEX COLORS
  if (is.null(V(net)$color)) {
    if (is.null(V(net)$group)) {
      V(net)$color <- rgb(0.5, 0.6, 1)
    } else {
      n <- length(unique(V(net)$group))
      if (n<=8) {
        col=brewer.pal(n, "Set2")
      } else {
        col=rainbow(n)
      }
      names(col) <- unique(V(net)$group)
      V(net)$color <- col[V(net)$group]
    }
  }
  V(net)$color <- sapply(V(net)$color, function(x) {add.alpha(x, alpha.v)})
  
  # VERTEX FRAME COLORS
  alpha.v <- min(alpha.v+0.3, 1)
  if (is.null(V(net)$frame.color)) {
    if (is.null(V(net)$sec.group)) {
      V(net)$frame.color <- V(net)$color
    } else {
      n <- length(unique(V(net)$sec.group))
      col=rainbow(n)
      names(col) <- unique(V(net)$sec.group)
      V(net)$frame.color <- col[V(net)$sec.group]
    }
  }
  V(net)$frame.color <- sapply(V(net)$color, function(x) {add.alpha(x, alpha.v)})
  
  # OTHER VERTEX ATTRIBUTES
  V(net)$label.color <- rgb(0, 0, .2, 1)
  V(net)$label.cex=2
  
  #EDGE ATTRIBUTES
  E(net)$color <- add.alpha(rgb(.5, .4, .5) , alpha.e)
  E(net)$lty <- ifelse(get.edge.attribute(net, name="value") > 0, 1, 3)
  E(net)$width <- scale.e(E(net)$value)
  
  ### PLOT
  if (is.null(layout.path)) {
    tkp = tkplot(net, layout <- layout.fruchterman.reingold, vertex.size = 5*scale.v(V(net)$weight), 
                 canvas.width=1000, canvas.height=1000)
    layout <- tkplot.getcoords(tkp)
    tk_off()
  } else {
    load(layout.path)
  }
  
  layout <- fit.to.layout(lay=layout)
  plot.igraph(net, layout=layout, vertex.size = scale.v(V(net)$weight), margin=0, rescale=F, xlim=c(0, 1), ylim=c(0, 1))
  
  # LEGEND
  if (!is.null(V(net)$group)) {
    legend(legend_x, legend_y, legend=unique(V(net)$group), pch=21, pt.bg = unique(V(net)$color), 
           pt.lwd = 0.1, bty="n", col = rgb(0.5, 0.5, 0.5, alpha = 0.3), y.intersp=1, cex=cex_legend)
  }
}

###########################33
#plot cooccur for one sample
###########################
plot.for.each.sample <- function(cooccurance, layout.path, weight, alpha.v=0.5, alpha.e=0.2, 
                                 scale.v = function(x) 10*x^0.5, scale.e = function(x) 20*(abs(x)^0.5), 
                                 legend_x= 0.7, legend_y=0.1, cex_legend=2.4) {
  net = copy(cooccurance)
  net$vertex$color = rgb(0.5, 0.6, 1)
  plot.graph(net, layout.path=layout.path, alpha.v=alpha.v-0.1, alpha.e=alpha.e, 
             scale.v=scale.v, scale.e=scale.e, legend_x=legend_x, legend_y=legend_y, cex_legend=cex_legend)
  net$vertex$weight = weight
  net$vertex$color = rgb(1, 0.4, 0.5)
  par(new=TRUE)
  plot.graph(net, layout.path=layout.path, alpha.v=alpha.v, alpha.e=0, 
             scale.v=scale.v, scale.e=scale.e, legend_x=legend_x, legend_y=legend_y, cex_legend=cex_legend)
}