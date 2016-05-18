
#--- loading case beta diversity ---
FamCaseM <- FamilyCase
ddHSN <- bcdist(FamCaseM)
FamCaseMMax <- FamCaseM [, which(colMaxs(FamCaseM) > 0)]
#--- end loading beta diversity ---

## create better names
rownames(Family) <- str_replace(rownames(Family), ".S.+", "_ctr")
rownames(Genus) <- str_replace(rownames(Genus), ".S.+", "_ctr")
rownames(Species) <- str_replace(rownames(Species), ".S.+", "_ctr")
rownames(Otu) <- str_replace(rownames(Otu), ".S.+", "_ctr")



#---Loading case and control (family, genus, species,otu, meta data, alpha diversity)---
Load <- function (FamInpCase, GenInpCase, SpeInpCase, OtuInpCase, MetaCaseCsv, AlphaDivInpCase, 
                  FamInpCtrl, GenInpCtrl, SpeInpCtrl, OtuInpCtrl, MetaCtrlCsv, AlphaDivInpCtrl)
{
  
  Family <- UniteMatrices(read_qiime_sum_feats (FamInpCase), read_qiime_sum_feats (FamInpCtrl))
  Genus <- UniteMatrices(read_qiime_sum_feats (GenInpCase), read_qiime_sum_feats (GenInpCtrl))
  Species <- UniteMatrices(read_qiime_sum_feats (SpeInpCase), read_qiime_sum_feats (SpeInpCtrl))
  Otu <- UniteMatrices(read_qiime_otu_table_no_tax (OtuInpCase), read_qiime_otu_table_no_tax (OtuInpCtrl))
  AlphaDiv <- UniteMatrices(LoadAlphaDiv(AlphaDivInpCase), LoadAlphaDiv(AlphaDivInpCtrl))
  
  # percent OTU
  Otup <- 100 * Otu / rowSums(Otu)
  
  ######### insert load alpha diversity as vector in list
  # load meta data
  MetaTable <- rbind(read.csv(MetaCaseCsv), read.csv(MetaCtrlCsv))
  
  list(Family=Family, Genus=Genus, Species=Species, Otu=Otu, Otup=Otup, Meta=MetaTable, AlphaDiv=AlphaDiv)
}
















