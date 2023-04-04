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


numCores <- 124
#registering clusters, can set a smaller number using numCores-1

registerDoParallel(numCores)

# Forever...

#load asymptotic d for two parameter distributions

#set the stop criterion
criterionset=1e-10

kurtWeibull<- read.csv(("kurtWeibull_28260.csv"))

allkurtWeibull<-unlist(kurtWeibull)

batchsize=100

samplesize=2048*2
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

#input the d value table previously generated
d_values<- read.csv(("d_SWA.csv"))
I_values<- read.csv(("I_SWA.csv"))
Imoments_values<- read.csv(("Imoments_SWA.csv"))
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
  n <- round(batchnumber*2048*9*3/100)
  (n%%10)==0
  # maximum order of moments
  morder <- 4
  #large sample size (approximating asymptotic)
  largesize<-round(batchnumber*2048*9/100)
  
  #generate quasirandom numbers based on the Sobol sequence
  quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                       mixed = FALSE, method = "C", start = 1)
  
  quasiuni<-quasiunisobol
  
  quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
  quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
  quasiuni_sorted4 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))
  
  
  
  samplesize=2048*2

  orderlist1_AB20<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=16,dimension=2)
  orderlist1_AB20<-orderlist1_AB20[1:largesize,]
  orderlist1_AB30<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=16,dimension=3)
  orderlist1_AB30<-orderlist1_AB30[1:largesize,]
  orderlist1_AB40<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=16,dimension=4)
  orderlist1_AB40<-orderlist1_AB40[1:largesize,]
  
  orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2,size=largesize,interval=16,dimension=2)
  orderlist1_AB2<-orderlist1_AB2[1:largesize,]
  orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3,size=largesize,interval=16,dimension=3)
  orderlist1_AB3<-orderlist1_AB3[1:largesize,]
  orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4,size=largesize,interval=16,dimension=4)
  orderlist1_AB4<-orderlist1_AB4[1:largesize,]
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    
    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1476],imoments1[1477],imoments1[1478])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[1345:1372])-targetm)/momentssd[1],
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[1396:1417])-targetvar)/momentssd[2],
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[1435:1448])-targettm)/momentssd[3],
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[1459:1468])-targetfm)/momentssd[4])
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[1345:1372]),
                 second=c(imoments1[2],momentsx[2],imoments1[1396:1417]),
                 third=c(imoments1[3],momentsx[3],imoments1[1435:1448]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[1459:1468])
    )
    medianmoments<-c(imoments1[1352],imoments1[1403],imoments1[1440],imoments1[1464])
    standardizedm<-c(imoments1[1352]/momentssd[1],imoments1[1403]/momentssd[2],imoments1[1440]/momentssd[3],imoments1[1464]/momentssd[4])
    all1<-(c(kurtx,skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }

  write.csv(SEbataches,paste("finite_Weibull_bootstrapsize_raw_SWA",batchnumber,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)

  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:39])^2))
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,40:63])^2))
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,64:79])^2))
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,80:91])^2))
  
  AB1_mean<-abs(colMeans((SEbataches[,7:39])))
  
  AB1_var<-abs(colMeans((SEbataches[,40:63])))
  
  AB1_tm<-abs(colMeans((SEbataches[,64:79])))
  
  AB1_fm<-abs(colMeans((SEbataches[,80:91])))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,108])
  
  samplemean_SE1<-samplemeansd_unscaled1
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,138])
  
  samplevar_SE1<-samplevarsd_unscaled1
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,162])
  
  sampletm_SE1<-sampletmsd_unscaled1
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,178])
  
  samplefm_SE1<-samplefmsd_unscaled1
  
  ratiosamplemean1<-c(SEbatachesmean[108])/SEbatachesmean[109]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,108])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1
  
  ratiosamplevar1<-c(SEbatachesmean[138])/SEbatachesmean[139]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,138])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1
  
  ratiosampletm1<-c(SEbatachesmean[162])/SEbatachesmean[163]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,162])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1
  
  ratiosamplefm1<-c(SEbatachesmean[178])/SEbatachesmean[179]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,178])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1
  
  ratiomean1<-c(SEbatachesmean[104:136])/SEbatachesmean[109]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,104:136]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,104:136])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1
  
  ratiovar1<-SEbatachesmean[137:160]/SEbatachesmean[139]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,137:160]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,137:160])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1
  ratiotm1<-SEbatachesmean[161:176]/SEbatachesmean[163]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,161:176]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,161:176])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1
  
  ratiofm1<-SEbatachesmean[177:188]/SEbatachesmean[179]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,177:188]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,177:188])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1
  
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

  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,108])
  
  samplemean_SE1<-samplemeansd_unscaled1
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,138])
  
  samplevar_SE1<-samplevarsd_unscaled1
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,162])
  
  sampletm_SE1<-sampletmsd_unscaled1
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,178])
  
  samplefm_SE1<-samplefmsd_unscaled1
  
  ratiosamplemean1<-c(SEbatachesmean[108])/SEbatachesmean[109]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,108])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1
  
  ratiosamplevar1<-c(SEbatachesmean[138])/SEbatachesmean[139]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,138])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1
  
  ratiosampletm1<-c(SEbatachesmean[162])/SEbatachesmean[163]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,162])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1
  
  ratiosamplefm1<-c(SEbatachesmean[178])/SEbatachesmean[179]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,178])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1
  
  ratiomean1<-c(SEbatachesmean[104:136])/SEbatachesmean[109]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,104:136]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,104:136])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1
  
  ratiovar1<-SEbatachesmean[137:160]/SEbatachesmean[139]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,137:160]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,137:160])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1
  ratiotm1<-SEbatachesmean[161:176]/SEbatachesmean[163]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,161:176]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,161:176])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1
  
  ratiofm1<-SEbatachesmean[177:188]/SEbatachesmean[179]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,177:188]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,177:188])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1
  
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
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatch_bias_Monte_SE,paste("finite_Weibull_bootstrapsize_raw_SWA_error",samplesize,".csv", sep = ","), row.names = FALSE)


