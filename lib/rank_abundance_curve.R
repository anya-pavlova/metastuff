library(reshape2)
library(plyr)

GetRankAbundenceCurveBasePlot <- function(top.by.sample.dt){
  RAC <- ggplot(top.by.sample.dt) + 
    geom_line(aes(x = rank, y = log(1 + abundance), group = sample), 
              alpha=0.03, size=1.5, color="black") +
    theme_bw() + xlab("rank") + ylab("abundance") +
    theme(text = element_text(size=6))
}

GetHighlitedRankPlot <- function(base.plot, top.by.sample.dt, sample.id){
  base.plot + geom_line(data = top.by.sample.dt[sample==sample.id],
                        aes(x = rank, y = log(1 + abundance)), color="red")
}

GetRankDataTable <- function(top.by.sample){
  rank.limit <- 20
  top.by.sample.dt <- as.data.table(ldply(.data = top.by.sample,
                                          .fun = function(x) unname(x)))
  setnames(top.by.sample.dt, "X1", "sample")
  
  top.by.sample.dt <- melt(top.by.sample.dt, id.vars = "sample")
  top.by.sample.dt$variable <- NULL
  top.by.sample.dt[,rank:=1:.N,by=sample]
  setnames(top.by.sample.dt, c("sample", "abundance","rank"))
  
  top.by.sample.dt <- top.by.sample.dt[,.SD[1:rank.limit],by=sample]
  top.by.sample.dt
}

