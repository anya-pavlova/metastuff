################################
# stat. tests
################################

THRESH_MAX_PVALUE <- 0.05     #0.1

# названия образцов интересующих нас для исследования и названия образцов контроля
#tags_ctrl <- sort(rownames(sp_ctrl))
#tags_vil <- sort(rownames(sp)[grep("^V", rownames(sp))])
#tags_Y89 <- sort(rownames(sp)[grep("^Y8|^Y9", rownames(sp))]) #for Y8 and Y9, yana and Borya 

################################
# families
################################

ff <- fam_all[c(tags_vil, tags_ctrl),] #матрица Только из вилюйцев и контроля


manes_farm <- sort(c(rownames(fam))) #вектор имен образцов
manes_farm_ctrl <- sort(c(rownames(fam_ctrl))) # вектор имен контроля

#colnames(ff) <- paste("o_", unlist(data.frame(strsplit(colnames(ff), "o_"))[2,]), sep="")

#ff <- ff[,which(colMaxs(ff) > 0)] # матрица из семейств без колонок в которых только 0
ff <- ff[,which(colMaxs(ff) > 0.1)] #  select top features
w <- wilcox_feats_both_dirs(ff, tags_vil, tags_ctrl, maxpv = THRESH_MAX_PVALUE, pairedt = FALSE)
if(length(w$greater) + length(w$less) > 0)
{
  print("Signif. diff. found.")
  wf <- convert_taxa_diff_to_table(w)
  wf2 <- convert_taxa_diff_to_table_DIV(w)  
  wf2[,"FDR adj. p-value"] <- round(wf2[,"FDR adj. p-value"], 3)
  write.table(wf2, file=paste("out/signif_fam.txt", sep=""), sep=" ", quote=F, row.names = F)    
} else {
  print("No signif. diff. found.")
}

# check a few
# one high in Vil.
b_vil <- ff[tags_vil, "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae"]
b_ctrl <- ff[tags_ctrl, "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae"]
boxplot(b_vil, b_ctrl)
sort(b_vil)
sort(fam["V9",])
# one high in controls
b_vil <- ff[tags_vil, "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae"]
b_ctrl <- ff[tags_ctrl, "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae"]
boxplot(b_vil, b_ctrl)



################################
# genera
################################
ff <- g_all[c(tags_vil, tags_ctrl),]
ff <- ff[,which(colMaxs(ff) > 0.1)]
colnames(ff) <- paste("f_", unlist(data.frame(strsplit(colnames(ff), "f_"))[2,]), sep="")
w <- wilcox_feats_both_dirs(ff, tags_vil, tags_ctrl, maxpv = THRESH_MAX_PVALUE, pairedt = FALSE)
if(length(w$greater) + length(w$less) > 0)
{
  print("Signif. diff. found.")
  wf <- convert_taxa_diff_to_table(w)
  wf2 <- convert_taxa_diff_to_table_DIV(w)  
  wf2[,"FDR adj. p-value"] <- round(wf2[,"FDR adj. p-value"], 3)
  write.table(wf2, file=paste("out/signif_g.txt", sep=""), sep=" ", quote=F, row.names = F)  
} else {
  print("No signif. diff. found.")
}


################################
# species
################################
ff <- sp_all[c(tags_vil, tags_ctrl),]
ff <- ff[,which(colMaxs(ff) > 0.1)]
colnames(ff) <- paste("o_", unlist(data.frame(strsplit(colnames(ff), "o_"))[2,]), sep="")
w <- wilcox_feats_both_dirs(ff, tags_vil, tags_ctrl, maxpv = THRESH_MAX_PVALUE, pairedt = FALSE)
if(length(w$greater) + length(w$less) > 0)
{
  print("Signif. diff. found.")
  wf <- convert_taxa_diff_to_table(w)
  wf2 <- convert_taxa_diff_to_table_DIV(w)  
  wf2[,"FDR adj. p-value"] <- round(wf2[,"FDR adj. p-value"], 3)
  write.table(wf2, file=paste("out/signif_org.txt", sep=""), sep=" ", quote=F, row.names = F)  
} else {
  print("No signif. diff. found.")
}


################################
# OTUs
################################
ff <- otup_all[c(tags_vil, tags_ctrl),]
ff <- ff[,which(colMaxs(ff) > 0.1)]
w <- wilcox_feats_both_dirs(ff, tags_vil, tags_ctrl, maxpv = THRESH_MAX_PVALUE, pairedt = FALSE)
if(length(w$greater) + length(w$less) > 0)
{
  print("Signif. diff. found.")
  w_bak <- w
  # add taxonomy to OTU names
  names(w$greater) <- paste(names(w$greater), taxotu[names(w$greater),], sep="_")
  names(w$less) <- paste(names(w$less), taxotu[names(w$less),], sep="_")
  # tabulated view with separating line
  wf2 <- convert_taxa_diff_to_table_DIV(w)  
  wf2[,"FDR adj. p-value"] <- round(wf2[,"FDR adj. p-value"], 3)
  # save separated diff-s to file
  write.table(wf2, file=paste("out/signif_otu.txt", sep=""), sep=" ", quote=F, row.names = F)
  # tabulated view without separator
  wf <- convert_taxa_diff_to_table(w_bak)  
} else {
  print("No signif. diff. found.")
}


