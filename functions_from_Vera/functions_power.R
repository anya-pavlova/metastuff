##################################
# functions for NeededSampleSize #
##################################

library(HMP)
library(futile.logger)

#adding features from another cohort and alphabetic ordering of features
merge.features <- function (group1,group2,eps){
  dataCase <- group1
  dataCntrl <- group2
  
  featCase <- colnames(dataCase)
  featCntrl <- colnames(dataCntrl)
  newfeatCntrl <- setdiff(featCase,featCntrl)
  newfeatCase <- setdiff(featCntrl,featCase)
  if(is.null(newfeatCase) && is.null(newfeatCntrl)) {
    dataCase <- dataCase[,order(colnames(dataCase))]
    dataCntrl <- dataCntrl[,order(colnames(DataControl))]
  }
  else{
    n<-nrow(dataCase)
    for (i in 1:length(newfeatCase)) {
      dataCase <- cbind(dataCase,rep(eps,n))
    }
    colnames(dataCase) <- c(featCase, newfeatCase)
    dataCase <- dataCase[,order(colnames(dataCase))]
    
    n <- nrow(dataCntrl)
    for (i in 1:length(newfeatCntrl)) {
      dataCntrl <- cbind(dataCntrl,rep(eps,n))
      i<-i+1
    } 
    colnames(dataCntrl) <- c(featCntrl, newfeatCntrl)
    dataCntrl <- dataCntrl[,order(colnames(dataCntrl))]
  }  
  list(dataCase=dataCase,dataCntrl=dataCntrl)
}

# Calculates the statistics, p-value and power of HMP comparisson of 2 groups,
# input: two matrices of counts with equal number of reads per sample, 
#        MC is the number of  Monte-Carlo experiments for HMP. In practice this should be at least 1,000.
# output: list containing staticstics value, p-value and power of MP comparisson of 2 groups
HMPcompare2gr <- function (gr1, gr2, MC = 1000){
  dataCase <- gr1
  dataCntrl <- gr2
  nrsCase <- rowSums(gr1)
  nrsCntrl <- rowSums(gr2)
  
  #comparing 2 groups with HMP
  mygroup <- list(dataCase, dataCntrl)
  compareCaseCntrl2 <- Xmcupo.sevsample(mygroup)
  print(compareCaseCntrl2)
  
  #calculating power
  fit.Case <- DM.MoM(dataCase)
  fit.Cntrl <- DM.MoM(dataCntrl)
  
  Nrs <- list(nrsCase, nrsCntrl)
  pi_2gr <- rbind(fit.Case$pi, fit.Cntrl$pi)
  group.theta <- rbind(fit.Case$theta,fit.Cntrl$theta)
  
  power <- MC.Xmcupo.statistics(Nrs, MC, fit.Cntrl$pi, pi_2gr, group.theta, "ha", 0.05)
  
  #summering results
  compareCaseCntrl2[3]  <- power
  names(compareCaseCntrl2)[3] <- "power"
  compareCaseCntrl2
}

#calculates the power of two groups comparison with HMP-package via different number for samples in one of them
#input: matrices of counts with equal number of reads per sample, step for increasing size of Case subsampling
#output: A matrix and  a plot of power vs Case subsampling size. The first column corresponding to numer of samples taken from the Case,
#the second column consists of corresponding value of power.
power.curve <- function (Case, Cntrl, step = 1, MC = 1500){
  # data <- merge.features(Case,Cntrl, 1e-10)
  group1 <- Case
  group2 <- Cntrl
  n1 <- nrow(group1)
  n2 <- nrow(group2)
  nc <- ncol(group1)
  s <- 3
  power <- c()
  nsamples <- c()
  i <- 1
  #takes n samples form each group and calculates power for them
  while (s <= n1)
  {
    usedrows1 <- c(1:s)
    usedrows2 <- c(1:n2)
    gr1 <- group1[usedrows1,]
    gr2 <- group2[usedrows2,]
    
    #checks wether there are unpresented in subgroups features
    for (l in 1:nc) if (sum(gr1[,l]) == 0) gr1[1,l] <- 1e-15
    for (l in 1:nc) if (sum(gr2[,l]) == 0) gr2[1,l] <- 1e-15
    
    #gr2[,which(colSums(gr2) == 0)] <- 1e-15
    
    
    power[i] <- HMPcompare2gr(gr1,gr2, MC)[[3]]
    nsamples[i] <- s 
    print (c(nsamples[i], power[i]))
    i <- i+1
    s <- s + step
  }
  if(s != n1 + step)
  {
    power[i] <- HMPcompare2gr(group1,group2, MC)[[3]]
    nsamples[i] <- n1
  }
  
  res <- matrix(c(nsamples, power), ncol = 2)
  colnames(res) <- c("number of samples", "power")
  plot(res[,1], res[,2])
  res
}


### Test functions

#calculates needed number of samples in group1 to compare it with group2 with some needed power, when all data is known.
samp_size <- function(group1, group2, needed_power = 0.8, MC = 1500){
  samp_pow <- power.curve(group1, group2, 1, MC)
  n_samp1 <- nrow(samp_pow)
  if (samp_pow[n_samp1,2] < needed_power) stop ("not enough samples to make a decision: power does not exeeds needed value")
  
  #flaggs the beginning of power curve's reliable part 
  flag <- n_samp1
  i <- n_samp1-1
  while (samp_pow[i,2] <= samp_pow[i+1,2] && i > 0){
    flag <- i+1
    i <- i-1
  }
  if (flag >= n_samp1-1) stop("not enough samples to make a desicion: unreliable power curve")
  
  #find the minimal sufficient number of samples
  i <- n_samp1
  while (samp_pow[i]>= needed_power && i>=flag) i <- i-1
  res <- samp_pow[i,1]
}
# test data for adding families 
#dataCase <- cbind(rep(1,5),rep(2,5),rep(3,5))
#dataCntrl <- cbind (rep(4,5),rep(5,5),rep(6,5),rep(7,5))
#colnames(dataCase) <- c("asdtsdg","bafgadfh","csdhs")
#colnames(dataCntrl) <- c("baadg","d","asdtsdg","e")

