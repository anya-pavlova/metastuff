
# draw PCoA biplot

dd <- bcdist(g)
res <- pcoa(dd)
# percentage of variance explained by each axis
perc_expl <- paste(round(100*(res$values$Relative_eig)[1:10], 1), "%", sep="")
pn <- paste("PC", 1:length(perc_expl), sep="")
perc_expl <- paste(pn, perc_expl, sep=": ")
perc_expl
tt <- g[,which(colMaxs(g) > 10)] # draw arrows for major only
#tt <- g

#pdf("graphs/biplot_123.pdf", width=15, height=7)
par(mfrow=c(1,2))
cg <- rownames(g)
#cg <- df_ctrl_case[,"tag",drop=F]
#cg[tags_ctrl,1] <- "*"
#cg[tags_crohn_stool,1] <- "S"
#cg[tags_crohn_il,1] <- "I"
#cg <- cg[,1]
length(cg)
nrow(tt)

biplot.pcoa.custom(res, tt, col=c("#0000ff", "black") , cex = c(0.8, 0.7), plot.axes=c(1,2), rn=cg, xlab_cus=perc_expl[1], ylab_cus=perc_expl[2])
#biplot.pcoa.custom(res, tt, col=c("#0000ff", "black"), cex = c(0.8, 0.7), plot.axes=c(1,3), rn=cg, xlab_cus=perc_expl[1], ylab_cus=perc_expl[3] )
#dev.off()
