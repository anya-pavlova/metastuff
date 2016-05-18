LoadPicrustBySample <- function(projects.directory.path, table.projects, project.ids, file.name="metagenome_predictions.tab"){
    predictions <- alply(unique(table.projects$folder), 1, function (x) {
      path <- file.path(projects.directory.path, x, file.name)
      samples <- table.projects[folder==x]$sample
      message(path)
      data.predictions <- fread(path, header=T, skip=1, sep = '\t', select=c('#OTU ID', samples))
    })
    #merge
    merged.predictions = predictions[[1]]
    for (prediction in predictions[-1]) {
      merged.predictions <- merge(merged.predictions, prediction, by = '#OTU ID', all=T)
    }
    colnames(merged.predictions)[1] <- 'kegg_id'
    duplicated.cols <- str_extract(colnames(merged.predictions), 'i\\..*')
    duplicated.cols <- duplicated.cols[!is.na(duplicated.cols)]
    if (length(duplicated.cols)>0) {
      merged.predictions <- merged.predictions[, str_replace(duplicated.cols, 'i\\.', ''):=NULL]
    }
    colnames(merged.predictions) <- str_replace(colnames(merged.predictions), 'i\\.', '')
    return(merged.predictions)
}

picrust_analysis_vitamins <- function(pred_data, quiime.data, vitamins, kegg_ko){
  kegg_id_vitamins <- merge(kegg_ko, vitamins[, c('kegg_ko', 'vitamin'), with=F], by='kegg_ko')
  pred.data.norm <- pred_data[kegg_id %in% kegg_id_vitamins$kegg_id, ]
  pred.data.norm <- merge(pred.data.norm, kegg_id_vitamins, by='kegg_id', all.x=T, all.y=F)
  pred.data.norm <- pred.data.norm[, c('vitamin', data.groups$sample_id), with=F]
  pred.data.norm <- melt(pred.data.norm[, lapply(.SD, sum), by=vitamin], id.vars=1,
                         value.name='value', variable.name='sample_id')
  
  plot_id <- GetVitamintBasePlot(pred.data.norm)
  quiime.data[, {
    p_id <- GetVitaminPlot(plot_id, pred.data.norm, before, col=rgb(1, 0.9, 0.1, alpha=1))
    p_id <- GetVitaminPlot(p_id, pred.data.norm, after, col=rgb(1, 0.1, 0.4, alpha=1))
    ggsave(file.path(dir, paste0(before, '_', after, '_vitamin.png')), device='png', width=10, height=5, plot = p_id)
    cat("Create graph for ", dir, "\n")
  }, by = 1:nrow(quiime.data)]
  
  pred.data.norm[, level:=calculate.level(value), by=vitamin]
  setorderv(pred.data.norm, c('sample_id', 'vitamin'))
  
  quiime.data[, {
    write.table(pred.data.norm[sample_id %in% c(before, after)], 
                file.path(dir, paste0(before, '_', after, '_level_vitamin.csv')),
                row.names = F, quote = F)
  }, by = 1:nrow(quiime.data)]
}

calculate.level <- function(value){
  thresh <- quantile(value, c(.33, .66)) 
  level <- rep('norm', length(value))
  level[value > thresh[2]] <- 'high'
  level[value < thresh[1]] <- 'low'
  return(level)
}

picrust_analysis_butyrate <- function(pred_data, quiime.data, butyrate){
  pred.data.norm <- pred_data[kegg_id %in% butyrate$kegg_id, ]
  pred.data.norm <- merge(pred.data.norm, butyrate, by='kegg_id', all.x=T, all.y=F)
  pred.data.norm <- pred.data.norm[, c('intermediate', data.groups$sample_id), with=F]
  pred.data.norm <- melt(pred.data.norm[, lapply(.SD, sum), by=intermediate], 
                         id.vars=1, value.name='value', variable.name='sample_id')
  
  plot_id <- GetButyrateBasePlot(pred.data.norm)
  quiime.data[, {
    p_id <- GetButyratePlot(plot_id, pred.data.norm, before, col=rgb(1, 0.9, 0.1, alpha=1))
    p_id <- GetButyratePlot(p_id, pred.data.norm, after, col=rgb(1, 0.1, 0.4, alpha=1))
    ggsave(file.path(dir, paste0(before, '_', after, '_butyrate.png')), device='png', width=10, height=5, plot = p_id)
    message("Create butyrate plot for ", dir, "\n")
  }, by = 1:nrow(quiime.data)]
  
  pred.data.norm[, level:=calculate.level(value), by=intermediate]
  setorderv(pred.data.norm, c('sample_id', 'intermediate'))
  
  quiime.data[, {
    write.table(pred.data.norm[sample_id %in% c(before, after)], 
                file.path(dir, paste0(before, '_', after, '_level_byterate.csv')),
                row.names = F, quote = F)
  }, by = 1:nrow(quiime.data)]
  
  pred_data_sum <- unique(pred.data.norm[, value:=sum(value), by = sample_id], by='sample_id')
  plot_sum_id <- GetButyrateSumBasePlot(pred_data_sum)
  quiime.data[, {
    p_sum_id <- GetButyrateSumPlot(plot_sum_id, pred.data.norm, before, col=rgb(1, 0.9, 0.1, alpha=1))
    p_sum_id <- GetButyrateSumPlot(p_sum_id, pred.data.norm, after, col=rgb(1, 0.1, 0.4, alpha=1))
    ggsave(file.path(dir, paste0(before, '_', after, '_sum_butyrate.png')), device='png', width=5, height=5, plot = p_sum_id)
    message("Create butyrate sum plot for ", dir, "\n")
  }, by = 1:nrow(quiime.data)]
}


GetVitamintBasePlot <- function(dt){
  PiCrust <- ggplot(dt, aes(x=vitamin, y=value)) + ylab("percent of genes") +
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) +
    theme_bw()
}

GetVitaminPlot <- function(base.plot, dt, sample.id, col=rgb(0.9, 0.8, 0.1, alpha=1)){
  PiCrust <- base.plot + geom_point(data = dt[sample_id==sample.id], size=3, color=col, pch=8)
}

GetButyrateBasePlot <- function(pred.data.norm){
  max.x <- max(pred.data.norm$intermediate)
  B.plot <- ggplot(pred.data.norm, aes(x=intermediate, y=value)) + geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) +
    ylab("percent of genes") + xlab("intermediate of butyrate metabolism")  + 
    scale_x_continuous(breaks= seq(1, max.x, by=1)) + theme_bw()
}

GetButyratePlot <- function(base.plot, pred.data.norm, id, col=rgb(0.9, 0.8, 0.1, alpha=1)){
  base.plot + geom_point(data = pred.data.norm[sample_id==id,], size=3, color=col, pch=8)
}

GetButyrateSumBasePlot <- function(pred.data.norm){
  B.plot <- ggplot(pred.data.norm, aes(x=1, y=value)) + ylab("percent of genes") + 
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) + theme_bw() + scale_x_continuous(breaks= seq(0,2,by=3))
}

GetButyrateSumPlot <- function(base.plot, pred.data.norm, id, col=rgb(0.9, 0.8, 0.1, alpha=1)){
  base.plot + geom_point(data = pred.data.norm[sample_id==id,], size=3, color=col, pch=8)
}
###################################

PiCrustVitaminGroup <- function(meta.1, meta.2,plot.name.1,plot.name.2, group.name ){
  
  pred.data.norm <- pred_table[kegg_id %in% kegg_id_vitamins$kegg_id, ]
  pred.data.norm <- merge(pred.data.norm, kegg_id_vitamins, by='kegg_id', all.x=T, all.y=F)
  pred.data.norm.1 <- pred.data.norm[, c('vitamin', meta.1), with=F]
  pred.data.norm.1 <- melt(pred.data.norm.1[, lapply(.SD, sum), by=vitamin], id.vars=1)
  
  pred.data.norm <- pred_table[kegg_id %in% kegg_id_vitamins$kegg_id, ]
  pred.data.norm <- merge(pred.data.norm, kegg_id_vitamins, by='kegg_id', all.x=T, all.y=F)
  pred.data.norm.2 <- pred.data.norm[, c('vitamin', meta.2), with=F]
  pred.data.norm.2 <- melt(pred.data.norm.2[, lapply(.SD, sum), by=vitamin], id.vars=1)
  
  PiCrust1 <- ggplot(pred.data.norm.1, aes(x=vitamin, y=value)) +
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) +
    theme_bw()+
    ggtitle(plot.name.1)
  
  PiCrust2 <- ggplot(pred.data.norm.2, aes(x=vitamin, y=value)) +
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) +
    theme_bw()+
    ggtitle(plot.name.2)
  
  SumPiCrust <- grid.arrange(PiCrust1, PiCrust2, ncol = 2)
  ggsave(paste0('/home/anna/metagenome/16s_study/out_ee_ka/PI-CRUST/', group.name, "_PI_CRUST",'.png'),width = 7, height= 5, SumPiCrust)
}

PiCrustButyrateGroup <- function(meta.data.1, meta.data.2, group.name, pred_data, butyrate,
                                 plot.name.1, plot.name.2){
  
  pred.data.norm <- pred_data[kegg_id %in% butyrate$kegg_id, ]
  pred.data.norm <- merge(pred.data.norm, butyrate, by='kegg_id', all.x=T, all.y=F)
  pred.data.norm.1 <- pred.data.norm[, c('intermediate', meta.data.1), with=F]
  pred.data.norm.1 <- melt(pred.data.norm.1[, lapply(.SD, sum), by=intermediate], id.vars=1, value.name='value', 
                           variable.name='sample_id')
  
  max.x.1 <- max(pred.data.norm.1$intermediate)
  B.plot.1 <- ggplot(pred.data.norm.1, aes(x=intermediate, y=value)) + 
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) +
    ylab("percent of genes") + xlab("intermediate of butyrate metabolism")  + 
    scale_x_continuous(breaks= seq(1, max.x.1, by=1)) + theme_bw()+
    ggtitle(plot.name.1)
  
  pred.data.norm <- pred_data[kegg_id %in% butyrate$kegg_id, ]
  pred.data.norm <- merge(pred.data.norm, butyrate, by='kegg_id', all.x=T, all.y=F)
  pred.data.norm.2 <- pred.data.norm[, c('intermediate', meta.data.2), with=F]
  pred.data.norm.2 <- melt(pred.data.norm.2[, lapply(.SD, sum), by=intermediate], id.vars=1, value.name='value', 
                           variable.name='sample_id')
  
  max.x.2 <- max(pred.data.norm.2$intermediate)
  B.plot.2 <- ggplot(pred.data.norm.2, aes(x=intermediate, y=value)) + 
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) +
    ylab("percent of genes") + xlab("intermediate of butyrate metabolism")  + 
    scale_x_continuous(breaks= seq(1, max.x.2, by=1)) + theme_bw()+
    ggtitle(plot.name.2)
  
  SumPiCrust <- grid.arrange(B.plot.1, B.plot.2, ncol = 2)
  ggsave(paste0('/home/anna/metagenome/16s_study/out_ee_ka/PI-CRUST/', group.name, "_butyrat",'.png'),width = 7, height= 5, SumPiCrust)
}

ButyratePlot <- function(pred.data.norm,outdir, group.name){
  max.x <- max(pred.data.norm$intermediate)
  B.plot <- ggplot(pred.data.norm, aes(x=intermediate, y=value)) + geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.1), pch=1) +
    ylab("percent of genes") + xlab("intermediate of butyrate metabolism")  + 
    scale_x_continuous(breaks= seq(1, max.x, by=1)) + theme_bw()
  ggsave(file.path(outdir, paste0(group.name,'_butyrate.png')), device='png', width=10, height=5, plot = plot_id)
}