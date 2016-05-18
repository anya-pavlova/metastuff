alphaDivCase <- (totalTable$AlphaDiv[which(rownames(totalTable$AlphaDiv) 
                                           %in% totalTable$Meta[which(totalTable$Meta[,"Type.1"] 
                                                                      %in% "case"),"samples_name"]),])

WriteTable(alphaDivCase, pathway$Case$OutdirCase, "AlphaDivCaseTbl")
AlphaDivMean <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], mean)
AlphaDivSd <- lapply(AlphaDivCaseL[match('AlphaDiversity', names(AlphaDivCaseL))], sd)
AlphaDivMeanSd <- c(AlphaDivMean, AlphaDivSd)
names(AlphaDivMeanSd) <- c("mean", "sd")

WriteTable (AlphaDivMeanSd, OutdirCase, "AlphaDivMeanAndSd") 

cairo_pdf((paste(OutdirCase, '/Graphs/alpha_boxplot.pdf', sep = "/")), width = 10,  height = 10)
boxplot(alphaDivCase$AlphaDiversity)
dev.off()

