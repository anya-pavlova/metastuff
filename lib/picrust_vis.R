GetPicrustBasePlot <- function(dt){
  PiCrust <- ggplot(dt, aes(x=vitamin, y=value)) +
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.2)) +
    theme_bw()
}

GetPicrusPlot <- function(base.plot, dt, sample.id){
  PiCrust <- base.plot + geom_point(data = dt[variable==sample.id], size=3, color=rgb(0.9, 0.6, 0.1, alpha=1))
}