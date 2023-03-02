#Copyright 2023 Tuban Lee

#These codes and manuscripts are under review in PNAS, please do not share them.

#If you are interested, please do not hesitate to contact me. Cooperation is also welcomed!

#require foreach and doparallel for parallel processing of bootstrap (not available for some types of computers)
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)
#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
if (!require("Rcpp")) install.packages("Rcpp")
library(Rcpp)
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)
if (!require("NRSReview")) install.packages("NRSReview_1.0.tar.gz", repos = NULL)
library(NRSReview)
if (!require("matrixStats")) install.packages("matrixStats")
library(matrixStats)


numCores <- detectCores()
#registering clusters, can set a smaller number using numCores-1

registerDoParallel(numCores)

# Forever...

#load asymptotic d for two parameter distributions

#set the stop criterion
criterionset=1e-30

kurtWeibull<- read.csv(("kurtWeibull_28260.csv"))

allkurtWeibull<-unlist(kurtWeibull)

batchsize=1000

samplesize=5400
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

#input the d value table previously generated
d_values<- read.csv(("d_SWA.csv"))
#Then, start the Monte Simulation

simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  library(randtoolbox)
  a=1

  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-c(kurtx=targetfm/(targetvar^(4/2)))
  skewx<-c(skewx=targettm/(targetvar^(3/2)))

  #bootsize for bootstrap approximation of the distributions of the kernal of U-statistics.
  n <- batchnumber*1.8*10^2
  (n%%10)==0
  # maximum order of moments
  morder <- 4
  #large sample size (approximating asymptotic)
  largesize<-1.8*10^4
  
  #generate quasirandom numbers based on the Sobol sequence
  quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                       mixed = FALSE, method = "C", start = 1)
  
  quasiuni<-quasiunisobol
  
  quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
  quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
  quasiuni_sorted4 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))
  
  
  
  samplesize=5400

  orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=8,dimension=2)
  orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=8,dimension=3)
  orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=8,dimension=4)
  
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    
    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    rqmoments1<-rqmoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,batch="auto",stepsize=100,criterion=criterionset)
    
    momentsx<-unbiasedmoments(x=sortedx)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),rqmoments1[153],rqmoments1[154],rqmoments1[155])

    allrawmoBias<-c(
      firstbias=(c(rqmoments1[seq(from=1, to=96, by=4)],rqmoments1[97:110])-targetm),
      secondbias=(c(rqmoments1[seq(from=2, to=96, by=4)],rqmoments1[111:124])-targetvar),
      thirdbias=(c(rqmoments1[seq(from=3, to=96, by=4)],rqmoments1[125:138])-targettm),
      fourbias=(c(rqmoments1[seq(from=4, to=96, by=4)],rqmoments1[139:152])-targetfm))

    all1<-t(c(kurtx,skewx,rqmoments1,targetall,momentsx,momentssd,allrawmoBias))

    SEbataches<-rbind(SEbataches,all1)
  }

  write.csv(SEbataches,paste("finite_Weibull_bootstrapsize_raw_SWA",batchnumber,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)

  SEbatachesmean <-apply(SEbataches, 2, calculate_column_mean)

  AB_mexp<-(which(abs(SEbatachesmean[seq(from=170, to=193, by=2)])==(min(abs(SEbatachesmean[seq(from=170, to=193, by=2)])))))[1]
  AB_mWeibull<-(which(abs(SEbatachesmean[seq(from=171, to=193, by=2)])==(min(abs(SEbatachesmean[seq(from=171, to=193, by=2)])))))[1]

  AB_varexp<-(which(abs(SEbatachesmean[seq(from=208, to=231, by=2)])==(min(abs(SEbatachesmean[seq(from=208, to=231, by=2)])))))[1]
  AB_varWeibull<-(which(abs(SEbatachesmean[seq(from=209, to=231, by=2)])==(min(abs(SEbatachesmean[seq(from=209, to=231, by=2)])))))[1]

  AB_tmexp<-(which(abs(SEbatachesmean[seq(from=246, to=269, by=2)])==(min(abs(SEbatachesmean[seq(from=246, to=269, by=2)])))))[1]
  AB_tmWeibull<-(which(abs(SEbatachesmean[seq(from=247, to=269, by=2)])==(min(abs(SEbatachesmean[seq(from=247, to=269, by=2)])))))[1]

  AB_fmexp<-(which(abs(SEbatachesmean[seq(from=284, to=307, by=2)])==(min(abs(SEbatachesmean[seq(from=284, to=307, by=2)])))))[1]
  AB_fmWeibull<-(which(abs(SEbatachesmean[seq(from=285, to=307, by=2)])==(min(abs(SEbatachesmean[seq(from=285, to=307, by=2)])))))[1]

  ABrank1<-c(AB_mexp=AB_mexp,AB_mWeibull=AB_mWeibull,AB_varexp=AB_varexp,AB_varWeibull=AB_varWeibull,AB_tmexp=AB_tmexp,AB_tmWeibull=AB_tmWeibull,AB_fmexp=AB_fmexp,AB_fmWeibull=AB_fmWeibull)

  ratiomean1<-c(SEbatachesmean[seq(from=3, to=98, by=4)],SEbatachesmean[99:112])/SEbatachesmean[99]

  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=3, to=98, by=4),99:112)]), 2, calculate_column_sd)

  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(seq(from=3, to=98, by=4),99:112)])/ratiomean1)

  meansd1<-apply(mean_SEbatachesmeanprocess, 2, calculate_column_sd)

  ratiovar1<-c(SEbatachesmean[seq(from=4, to=98, by=4)],SEbatachesmean[113:126])/SEbatachesmean[113]

  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=4, to=98, by=4),113:126)]), 2, calculate_column_sd)

  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=4, to=98, by=4),113:126)])/ratiovar1))

  varsd1<-apply(var_SEbatachesvarprocess, 2, calculate_column_sd)

  ratiotm1<-c(SEbatachesmean[seq(from=5, to=98, by=4)],SEbatachesmean[127:140])/SEbatachesmean[127]

  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=5, to=98, by=4),127:140)]), 2, calculate_column_sd)

  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=5, to=98, by=4),127:140)])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, calculate_column_sd)

  ratiofm1<-SEbatachesmean[ c(seq(from=6, to=98, by=4),141:154)]/SEbatachesmean[141]
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=6, to=98, by=4),141:154)]), 2, calculate_column_sd)

  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=6, to=98, by=4),141:154)])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, calculate_column_sd)

  allSE_unstan<-c(meansd_unscaled1=meansd_unscaled1,varsd_unscaled1=varsd_unscaled1,tmsd_unscaled1=tmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE_unstand<-c(meansd1=meansd1,varsd1=varsd1,tmsd1=tmsd1,fmsd1=fmsd1
  )

  mexp<-which(meansd1[seq(from=1, to=24, by=2)]==(min(meansd1[seq(from=1, to=24, by=2)])))[1]
  mWeibull<-which(meansd1[seq(from=2, to=24, by=2)]==(min(meansd1[seq(from=2, to=24, by=2)])))[1]
  varexp<-which(varsd1[seq(from=1, to=24, by=2)]==(min(varsd1[seq(from=1, to=24, by=2)])))[1]
  varWeibull<-which(varsd1[seq(from=2, to=24, by=2)]==(min(varsd1[seq(from=2, to=24, by=2)])))[1]
  tmexp<-which(tmsd1[seq(from=1, to=24, by=2)]==(min(tmsd1[seq(from=1, to=24, by=2)])))[1]
  tmWeibull<-which(tmsd1[seq(from=2, to=24, by=2)]==(min(tmsd1[seq(from=2, to=24, by=2)])))[1]
  fmexp<-which(fmsd1[seq(from=1, to=24, by=2)]==(min(fmsd1[seq(from=1, to=24, by=2)])))[1]
  fmWeibull<-which(fmsd1[seq(from=2, to=24, by=2)]==(min(fmsd1[seq(from=2, to=24, by=2)])))[1]

  rank1<-c(mexp=mexp,mWeibull=mWeibull,varexp=varexp,varWeibull=varWeibull,tmexp=tmexp,tmWeibull=tmWeibull,fmexp=fmexp,fmWeibull=fmWeibull)

  allresultsSE<-c(samplesize=samplesize,type=1,SEbatachesmean,allSE_unstan,allSSE_unstand,ABrank1,rank1)
}

write.csv(simulatedbatch_bias_Monte,paste("finite_Weibull_bootstrapsize_raw_SWA",samplesize,".csv", sep = ","), row.names = FALSE)


simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)

  a=1

  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<- read.csv(paste("finite_Weibull_bootstrapsize_raw_SWA",batchnumber,round(kurtx,digits = 1),".csv", sep = ","))

  SEbatachesmean <-apply(SEbataches, 2, calculate_column_mean)

  ratiomean1<-c(SEbatachesmean[seq(from=3, to=98, by=4)],SEbatachesmean[99:112])/SEbatachesmean[99]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=3, to=98, by=4),99:112)]), 2, se_sd)
  
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(seq(from=3, to=98, by=4),99:112)])/ratiomean1)
  
  meansd1<-apply(mean_SEbatachesmeanprocess, 2, se_sd)
  
  ratiovar1<-c(SEbatachesmean[seq(from=4, to=98, by=4)],SEbatachesmean[113:126])/SEbatachesmean[113]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=4, to=98, by=4),113:126)]), 2, se_sd)
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=4, to=98, by=4),113:126)])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  ratiotm1<-c(SEbatachesmean[seq(from=5, to=98, by=4)],SEbatachesmean[127:140])/SEbatachesmean[127]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=5, to=98, by=4),127:140)]), 2, se_sd)
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=5, to=98, by=4),127:140)])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  
  ratiofm1<-SEbatachesmean[ c(seq(from=6, to=98, by=4),141:154)]/SEbatachesmean[141]
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=6, to=98, by=4),141:154)]), 2, se_sd)
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=6, to=98, by=4),141:154)])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  
  allSE_unstandardized_se_sd<-c(meansd_unscaled1=meansd_unscaled1,varsd_unscaled1=varsd_unscaled1,tmsd_unscaled1=tmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE_unstandardized_se_sd<-c(meansd1=meansd1,varsd1=varsd1,tmsd1=tmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allresultsSE<-c(samplesize=samplesize,type=1,SEbatachesmean,se_mean_all1=se_mean_all1,allSE_unstandardized_se_sd=allSE_unstandardized_se_sd,allSSE_unstandardized_se_sd=allSSE_unstandardized_se_sd)

  allresultsSE
}

write.csv(simulatedbatch_bias_Monte_SE,paste("finite_Weibull_bootstrapsize_raw_SWA_error",samplesize,".csv", sep = ","), row.names = FALSE)


