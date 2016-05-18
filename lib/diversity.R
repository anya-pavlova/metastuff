
library(data.table)
library(plyr)

ReadAlpha <- function(path)
{
  data.alpha <- data.table(fread(path, colClasses = 'character', header=F))
  # renaming some columns
  setnames(data.alpha, "V1", "name")
  setnames(data.alpha, "V2", "seq_per_sample")
  setnames(data.alpha, "V3", "iteration")
  # 1st line is header
  data.alpha.header <- as.matrix(data.alpha[1]) 
  # everything below header are values
  data.alpha.values <- data.alpha[2:nrow(data.alpha)] 
  # convert "sequence per sample" to numbers because we load it as character
  data.alpha.values$seq_per_sample <- as.numeric(data.alpha.values$seq_per_sample) 
  # find max sequence per sample and take only those rows where it max
  max.seq <- max(data.alpha.values$seq_per_sample)
  names.to.keep <- data.alpha.values[seq_per_sample==max.seq]$name
  data.alpha.values.sub <- data.alpha.values[name %in% names.to.keep]
  
  # put into `numbers.col.names` names of columns with diversity values
  numbers.col.names <- names(data.alpha.values.sub)[4:ncol(data.alpha.values.sub)]
  # convert each column in `data.alpha.values.sub` to numeric
  data.alpha.values.sub.numbers <- data.alpha.values.sub[, lapply(.SD, as.numeric), .SDcols=numbers.col.names]
  # take mean of each column of `data.alpha.values.sub.numbers`
  data.alpha.values.sub.numbers.mean <- data.alpha.values.sub.numbers[, lapply(.SD, mean), .SDcols=numbers.col.names]
  
  # now check that means are correct
  # V4 variable always exists, it is the first diversity value
  if (mean(as.numeric(data.alpha.values.sub$V4)) != data.alpha.values.sub.numbers.mean$V4){
    warning("means for alpha diversity are wrong calculated. Please debug")
    return("error")
  }
  
  # take sample names from header
  data.alpha.header.samples <- data.alpha.header[4:length(data.alpha.header)]
  # now transpose means and cbind samples names column. Order is preserved
  data.alpha.values.sub.numbers.mean.t <- t(data.alpha.values.sub.numbers.mean)
  data.diversity <- data.table(cbind(data.alpha.values.sub.numbers.mean.t, data.alpha.header.samples))
  setnames(data.diversity, c("diversity", "sample"))
  data.diversity$diversity <- as.numeric(data.diversity$diversity)
  data.diversity
}

ChooseAlphaSamples <- function(alpha.div.project, data.groups)
{
  alpha.div.project.proj <- alpha.div.project[sample %in% data.groups$sample_id,]
}  

GetAlphaDivBasePlot <- function(data){
  require(ggplot2)
  base.plot <- ggplot(data, aes(x=group, y=diversity_scaled)) + ylab("diversity")+
    xlab("")+ 
    geom_point(size=3, color=rgb(0.2, 0.7, 0.8, alpha=0.5), pch=1) + 
    theme_bw() + scale_y_continuous(breaks=seq(1, max(data$diversity_scaled), by=1))
}

GetAlphaDivSamplePlot <- function(baseplot, data, sample.name){
  require(ggplo2)
  setkey(alpha.div.project, 'sample')
  baseplot + 
    geom_point(data = data[sample.name,], size=3, pch=8, col=rgb(1, 0.1, 0.4, alpha=1))
}

MakeAlphaBoxPlot <- function(alpha.data, plotname){
  require(ggplot2)
  ggplot(alpha.data, aes(x=group, y=diversity, group=group, colour= factor(group))) + 
    geom_boxplot() +
    theme_bw()+
    scale_x_discrete(limit = c("After", "Before"))+
    guides(color=guide_legend(title="after/before"))+
    ggtitle(plotname)
  ggsave(paste0('out_ee_ka/alpha_div/', plotname,'.png'),  width = 7, height= 5)
}
# alpha.div.project.g1
# alpha.div.project.g2
# alpha.div.project 
# 
# setkey(alpha.div.project.g1, 'sample')
# setkey(alpha.div.project.g2, 'sample')
# # 
# # AlpaDivRes("01-009-02","06-055-03")
# # AlpaDivRes <- function(id.semple.g1, id.semple.g2)
# # {
#    id.semple.g1="1.121.1"
#    id.semple.g2="1.121.2"
#    setkey(alpha.div.project.g1, "sample")
#    setkey(alpha.div.project.g2, "sample")
#    ind.sempl.af <-alpha.div.project.g1[id.semple.g1]
#    ind.sempl.bef <-alpha.div.project.g2[id.semple.g2]
#    cairo_pdf((paste('/home/anna/metagenome/16s_study/figures/alpha_boxplotTEST.pdf',
#                     sep = "/")), width = 10,  height = 10)
#    
#    alpha.div.project.mean.g1 <- mean(alpha.div.project.g1$diversity)
# #   alpha.div.project.mean.g2 <- mean(alpha.div.project.g2$diversity)
# #   
#    alpha.div.project.sd.g1 <- sd(alpha.div.project.g1$diversity)
#    alpha.div.project.sd.g2 <- sd(alpha.div.project.g2$diversity)
# 
#   ggplot(alpha.div.project.g1, aes(x=group, y=diversity, group=group)) + 
#     geom_boxplot() +
#     theme_bw()+
#     scale_x_discrete(limit = c("After", "Before"))+
#     #theme(panel.border = element_blank())+  
#     geom_point(data=ind.sempl.af, size = 5, color='darkorchid4', alpha=0.6)+
#     geom_text(data = ind.sempl.af, col="black", label="you before")+
#     # geom_point(data=ind.sempl.bef, size = 5, color='darkorchid4', alpha=0.6)+
#     # geom_text(data = ind.sempl.bef, col="black", label="you after")+
#     #stat_summary(fun.data=mean_sdl, mult=1, 
#     #            geom="errorbar", color="red", width=0.2) +
#     geom_errorbar(aes(ymin=diversity-alpha.div.project.sd.g1, ymax=diversity+alpha.div.project.sd.g1), width=.2,
#                   position=position_dodge(.9))   
#   
#     ggsave('/home/anna/metagenome/16s_study/figures/alpha_boxplotTEST.pdf',  width = 7, height= 5)
#      #geom_errorbar(aes(ymin=mean-alpha.div.project.sd.g1,ymax=mean+sd),linetype = 3,width = 0.25)
#      #stat_boxplot(geom ='errorbar')
#      #facet_grid(. ~ group) #разнести по двум графикам
# 
# # }

