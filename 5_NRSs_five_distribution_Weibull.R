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


#bootsize for bootstrap approximation of the distributions of the kernel of U-statistics.
n <- 1.8*10^4
(n%%10)==0
# maximum order of moments
morder <- 4

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)
quasiuni<-rbind(quasiunisobol)

quasiunisobol<-c()

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...

asymptotic_n <- 1.8*10^6
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (asymptotic bias)
largesize<-1.8*10^6

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_asymptotic<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

quasiuni_sorted2_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4_asymptotic <- na.omit(rowSort(quasiuni_asymptotic, descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...

quasiuni_asymptotic<-Sort(quasiuni_asymptotic[,1])

orderlist1_AB2_asymptotic<-createorderlist(quni1=quasiuni_sorted2_asymptotic,size=largesize,interval=8,dimension=2)
orderlist1_AB3_asymptotic<-createorderlist(quni1=quasiuni_sorted3_asymptotic,size=largesize,interval=8,dimension=3)
orderlist1_AB4_asymptotic<-createorderlist(quni1=quasiuni_sorted4_asymptotic,size=largesize,interval=8,dimension=4)

quasiuni_sorted2_asymptotic<-c()
quasiuni_sorted3_asymptotic<-c()
quasiuni_sorted4_asymptotic<-c()

#load asymptotic d for two parameter distributions
d_values<- read.csv(("d_SWA.csv"))
I_values<- read.csv(("I_SWA.csv"))

#set the stop criterion
criterionset=1e-30

kurtWeibull<- read.csv(("kurtWeibull_31150.csv"))
allkurtWeibull<-unlist(kurtWeibull)

simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtWeibull)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=allkurtWeibull[batchnumber]
  x<-c(dsWeibull(uni=quasiuni_asymptotic, shape=a/1, scale = 1))
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2_asymptotic)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  
  momentsx<-unbiasedmoments(x=sortedx)

  sortedx<-c()

  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])

  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm)/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm)/momentssd[4])
  
  medianmoments<-c(imoments1[104],imoments1[118],imoments1[132],imoments1[146])
  standardizedm<-c(imoments1[104]/momentssd[1],imoments1[118]/momentssd[2],imoments1[132]/momentssd[3],imoments1[146]/momentssd[4])
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_Weibull_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtWeibull),1],simulatedbatch_asymptoticbias[1:length(allkurtWeibull),c(193:ncol(simulatedbatch_asymptoticbias))]),paste("asymptotic_Weibull",largesize,".csv", sep = ","), row.names = FALSE)

samplesize=5400
batchsizebase=1000
orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=8,dimension=2)
orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=8,dimension=3)
orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=8,dimension=4)


batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)

  a=allkurtWeibull[batchnumber]

  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<-c()
  for (batch1 in c(1:batchsize)){

    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx)
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
    HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,batch="auto",stepsize=100,criterion=criterionset,boot=TRUE)
    
    momentsx<-unbiasedmoments(x=sortedx)

    sortedx<-c()

    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])

    allrawmoBias<-c(
      firstbias=(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm),
      secondbias=(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar),
      thirdbias=(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm),
      fourbias=(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("finite_Weibull_raw_SWA",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,193:246])^2))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,247:289])^2))/simulatedbatch_asymptoticbias[batchnumber,377]

  RMSE1_tm<-sqrt(colMeans((SEbataches[,290:332])^2))/simulatedbatch_asymptoticbias[batchnumber,378]

  RMSE1_fm<-sqrt(colMeans((SEbataches[,333:375])^2))/simulatedbatch_asymptoticbias[batchnumber,379]

  AB1_mean<-abs(colMeans((SEbataches[,193:246])))/simulatedbatch_asymptoticbias[batchnumber,376]

  AB1_var<-abs(colMeans((SEbataches[,247:289])))/simulatedbatch_asymptoticbias[batchnumber,377]

  AB1_tm<-abs(colMeans((SEbataches[,290:332])))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  AB1_fm<-abs(colMeans((SEbataches[,333:375])))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,189])

  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]

  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,190])

  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]

  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,191])

  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]

  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,192])

  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)

  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]

  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]

  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)

  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]

  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]

  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)

  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]

  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)

  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]

  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, unbiasedsd)

  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)

  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]

  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]

  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, unbiasedsd)

  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]

  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))

  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)

  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]

  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, unbiasedsd)

  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]

  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]

  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]

  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, unbiasedsd)

  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]

  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )

  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)


  allErrors
}

write.csv(simulatedbatch_ABSE,paste("Weibull_ABSSE.csv", sep = ","), row.names = FALSE)



simulatedbatch_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtWeibull[batchnumber]

  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<- read.csv(paste("finite_Weibull_raw_SWA",samplesize,round(kurtx,digits = 1),".csv", sep = ","))


  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)

  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("Weibull_ABSSE_error.csv", sep = ","), row.names = FALSE)





kurtgamma<- read.csv(("kurtgamma_31150.csv"))
allkurtgamma<-unlist(kurtgamma)

simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtgamma)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=allkurtgamma[batchnumber]
  x<-c(dsgamma(uni=quasiuni_asymptotic, shape=a/1, rate  = 1))
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2_asymptotic)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm)/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm)/momentssd[4])
  
  medianmoments<-c(imoments1[104],imoments1[118],imoments1[132],imoments1[146])
  standardizedm<-c(imoments1[104]/momentssd[1],imoments1[118]/momentssd[2],imoments1[132]/momentssd[3],imoments1[146]/momentssd[4])
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_gamma_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtgamma),1],simulatedbatch_asymptoticbias[1:length(allkurtgamma),c(193:ncol(simulatedbatch_asymptoticbias))]),paste("asymptotic_gamma",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtgamma))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtgamma[batchnumber]

  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<-c()
  for (batch1 in c(1:batchsize)){

    x<-c(dsgamma(uni=unibatch[,batch1], shape=a/1, rate  = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx)
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
    HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,batch="auto",stepsize=100,criterion=criterionset,boot=TRUE)
    
    momentsx<-unbiasedmoments(x=sortedx)
    
    sortedx<-c()
    
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
    
    allrawmoBias<-c(
      firstbias=(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm),
      secondbias=(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar),
      thirdbias=(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm),
      fourbias=(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
    
    SEbataches<-rbind(SEbataches,all1)
  }

  write.csv(SEbataches,paste("gamma_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)

  RMSE1_mean<-sqrt(colMeans((SEbataches[,193:246])^2))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,247:289])^2))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,290:332])^2))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,333:375])^2))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  AB1_mean<-abs(colMeans((SEbataches[,193:246])))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  AB1_var<-abs(colMeans((SEbataches[,247:289])))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  AB1_tm<-abs(colMeans((SEbataches[,290:332])))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  AB1_fm<-abs(colMeans((SEbataches[,333:375])))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_ABSE,paste("gamma_ABSSE.csv", sep = ","), row.names = FALSE)



simulatedbatch_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtgamma))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtgamma[batchnumber]

  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<- read.csv(paste("gamma_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))


  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("gamma_ABSSE_error.csv", sep = ","), row.names = FALSE)




kurtPareto<- read.csv(("kurtPareto_91210.csv"))
allkurtPareto<-unlist(kurtPareto)

simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtPareto)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=allkurtPareto[batchnumber]
  x<-c(dsPareto(uni=quasiuni_asymptotic, shape=a/1, scale = 1))

  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2_asymptotic)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm)/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm)/momentssd[4])
  
  medianmoments<-c(imoments1[104],imoments1[118],imoments1[132],imoments1[146])
  standardizedm<-c(imoments1[104]/momentssd[1],imoments1[118]/momentssd[2],imoments1[132]/momentssd[3],imoments1[146]/momentssd[4])
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_Pareto_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtPareto),1],simulatedbatch_asymptoticbias[1:length(allkurtPareto),193:ncol(simulatedbatch_asymptoticbias)]),paste("asymptotic_Pareto",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtPareto))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtPareto[batchnumber]

  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<-c()
  for (batch1 in c(1:batchsize)){

    x<-c(dsPareto(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx)
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
    HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,batch="auto",stepsize=100,criterion=criterionset,boot=TRUE)
    
    momentsx<-unbiasedmoments(x=sortedx)
    
    sortedx<-c()
    
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
    
    allrawmoBias<-c(
      firstbias=(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm),
      secondbias=(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar),
      thirdbias=(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm),
      fourbias=(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
    
    SEbataches<-rbind(SEbataches,all1)
  }

  write.csv(SEbataches,paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)

  RMSE1_mean<-sqrt(colMeans((SEbataches[,193:246])^2))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,247:289])^2))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,290:332])^2))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,333:375])^2))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  AB1_mean<-abs(colMeans((SEbataches[,193:246])))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  AB1_var<-abs(colMeans((SEbataches[,247:289])))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  AB1_tm<-abs(colMeans((SEbataches[,290:332])))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  AB1_fm<-abs(colMeans((SEbataches[,333:375])))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_ABSE,paste("Pareto_ABSSE.csv", sep = ","), row.names = FALSE)



simulatedbatch_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtPareto))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtPareto[batchnumber]

  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<- read.csv(paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))


  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("Pareto_ABSSE_error.csv", sep = ","), row.names = FALSE)




kurtlognorm<- read.csv(("kurtlognorm_31150.csv"))
allkurtlognorm<-unlist(kurtlognorm)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtlognorm)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=allkurtlognorm[batchnumber]
  x<-c(dslnorm(uni=quasiuni_asymptotic, meanlog =0, sdlog  = a/1))
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2_asymptotic)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm)/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm)/momentssd[4])
  
  medianmoments<-c(imoments1[104],imoments1[118],imoments1[132],imoments1[146])
  standardizedm<-c(imoments1[104]/momentssd[1],imoments1[118]/momentssd[2],imoments1[132]/momentssd[3],imoments1[146]/momentssd[4])
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_lognorm_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtlognorm),1],simulatedbatch_asymptoticbias[1:length(allkurtlognorm),193:ncol(simulatedbatch_asymptoticbias)]),paste("asymptotic_lognorm",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtlognorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)

  a=allkurtlognorm[batchnumber]

  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<-c()
  for (batch1 in c(1:batchsize)){

    x<-c(dslnorm(uni=unibatch[,batch1],  meanlog =0, sdlog = a/1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx)
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
    HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,batch="auto",stepsize=100,criterion=criterionset,boot=TRUE)
    
    momentsx<-unbiasedmoments(x=sortedx)
    
    sortedx<-c()
    
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
    
    allrawmoBias<-c(
      firstbias=(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm),
      secondbias=(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar),
      thirdbias=(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm),
      fourbias=(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
    
    SEbataches<-rbind(SEbataches,all1)
  }

  write.csv(SEbataches,paste("lognorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)

  RMSE1_mean<-sqrt(colMeans((SEbataches[,193:246])^2))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,247:289])^2))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,290:332])^2))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,333:375])^2))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  AB1_mean<-abs(colMeans((SEbataches[,193:246])))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  AB1_var<-abs(colMeans((SEbataches[,247:289])))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  AB1_tm<-abs(colMeans((SEbataches[,290:332])))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  AB1_fm<-abs(colMeans((SEbataches[,333:375])))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_ABSE,paste("lognorm_ABSSE.csv", sep = ","), row.names = FALSE)



simulatedbatch_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtlognorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtlognorm[batchnumber]

  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<- read.csv(paste("lognorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))


  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("lognorm_ABSSE_error.csv", sep = ","), row.names = FALSE)





kurtgnorm<- read.csv(("kurtgnorm_31150.csv"))
allkurtgnorm<-unlist(kurtgnorm)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtgnorm)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=allkurtgnorm[batchnumber]
  x<-c(dsgnorm(uni=quasiuni_asymptotic, shape=a/1, scale = 1))

  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2_asymptotic)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm)/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm)/momentssd[4])
  
  medianmoments<-c(imoments1[104],imoments1[118],imoments1[132],imoments1[146])
  standardizedm<-c(imoments1[104]/momentssd[1],imoments1[118]/momentssd[2],imoments1[132]/momentssd[3],imoments1[146]/momentssd[4])
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_gnorm_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtgnorm),1],simulatedbatch_asymptoticbias[1:length(allkurtgnorm),193:ncol(simulatedbatch_asymptoticbias)]),paste("asymptotic_gnorm",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtgnorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtgnorm[batchnumber]

  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<-c()
  for (batch1 in c(1:batchsize)){

    x<-c(dsgnorm(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx)
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
    HL1<-Hodges_Lehmann(x=sortedx,orderlist1_sorted2 = orderlist1_AB2)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,batch="auto",stepsize=100,criterion=criterionset,boot=TRUE)
    
    momentsx<-unbiasedmoments(x=sortedx)
    
    sortedx<-c()
    
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
    
    allrawmoBias<-c(
      firstbias=(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-targetm),
      secondbias=(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-targetvar),
      thirdbias=(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-targettm),
      fourbias=(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
    
    SEbataches<-rbind(SEbataches,all1)
  }

  write.csv(SEbataches,paste("gnorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)

  RMSE1_mean<-sqrt(colMeans((SEbataches[,193:246])^2))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,247:289])^2))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,290:332])^2))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,333:375])^2))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  AB1_mean<-abs(colMeans((SEbataches[,193:246])))/simulatedbatch_asymptoticbias[batchnumber,376]
  
  AB1_var<-abs(colMeans((SEbataches[,247:289])))/simulatedbatch_asymptoticbias[batchnumber,377]
  
  AB1_tm<-abs(colMeans((SEbataches[,290:332])))/simulatedbatch_asymptoticbias[batchnumber,378]
  
  AB1_fm<-abs(colMeans((SEbataches[,333:375])))/simulatedbatch_asymptoticbias[batchnumber,379]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_ABSE,paste("gnorm_ABSSE.csv", sep = ","), row.names = FALSE)



simulatedbatch_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtgnorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)


  a=allkurtgnorm[batchnumber]

  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  SEbataches<- read.csv(paste("gnorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))


  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,189])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,190])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,191])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,192])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiosamplemean1<-c(SEbatachesmean[189])/SEbatachesmean[110]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,189])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiosamplevar1<-c(SEbatachesmean[190])/SEbatachesmean[124]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,190])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  ratiosampletm1<-c(SEbatachesmean[191])/SEbatachesmean[138]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,191])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiosamplefm1<-c(SEbatachesmean[192])/SEbatachesmean[152]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,192])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  ratiomean1<-c(SEbatachesmean[c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/SEbatachesmean[110]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,376]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(3:13,seq(from=14, to=109, by=4),110:123,seq(from=169, to=184, by=4))])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatch_asymptoticbias[batchnumber,376]
  
  ratiovar1<-SEbatachesmean[c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]/SEbatachesmean[124]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,377]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=15, to=109, by=4),124:137,seq(from=170, to=184, by=4))])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,377]
  ratiotm1<-SEbatachesmean[c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]/SEbatachesmean[138]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=16, to=109, by=4),138:151,seq(from=171, to=184, by=4))])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,378]
  
  ratiofm1<-SEbatachesmean[c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]/SEbatachesmean[152]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,c(seq(from=17, to=109, by=4),152:165,seq(from=172, to=184, by=4))])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,379]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("gnorm_ABSSE_error.csv", sep = ","), row.names = FALSE)






