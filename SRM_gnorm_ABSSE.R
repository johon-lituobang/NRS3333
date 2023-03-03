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
morder <- 5

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)
quasiuni<-rbind(quasiunisobol)

quasiunisobol<-c()

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni[,1:4], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted5 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...

asymptotic_n <- 1.8*10^6
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 5
#large sample size (asymptotic bias)
largesize<-1.8*10^6

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_asymptotic<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

quasiuni_sorted2_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:4], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted5_asymptotic <- na.omit(rowSort(quasiuni_asymptotic, descend = FALSE, stable = FALSE, parallel = TRUE))

# Forever...

quasiuni_asymptotic<-Sort(quasiuni_asymptotic[,1])

orderlist1_AB2_asymptotic<-createorderlist(quni1=quasiuni_sorted2_asymptotic,size=largesize,interval=8,dimension=2)
orderlist1_AB3_asymptotic<-createorderlist(quni1=quasiuni_sorted3_asymptotic,size=largesize,interval=8,dimension=3)
orderlist1_AB4_asymptotic<-createorderlist(quni1=quasiuni_sorted4_asymptotic,size=largesize,interval=8,dimension=4)
orderlist1_AB5_asymptotic<-createorderlist(quni1=quasiuni_sorted5_asymptotic,size=largesize,interval=8,dimension=5)

quasiuni_sorted2_asymptotic<-c()
quasiuni_sorted3_asymptotic<-c()
quasiuni_sorted4_asymptotic<-c()
quasiuni_sorted5_asymptotic<-c()

#set the stop criterion
criterionset=1e-30


median_of_means <- function(x,korder=2) {
  n_blocks=ceiling(length(x)/korder)
  if (n_blocks > length(x)) {
    n_blocks <- ceiling(length(x) / 2)
  }
  indic <- sample(rep(1:n_blocks, each =korder),length(x), replace = FALSE)
  
  meangroup1<-split(x, indic)
  meangroup1<-as.data.frame(meangroup1)
  groupmean<-colMeans(meangroup1)
  return(median(groupmean))
}


kurtgnorm<- read.csv(("kurtgnorm_31150.csv"))
allkurtgnorm<-unlist(kurtgnorm)


samplesize=5400
batchsizebase=1000
orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=8,dimension=2)
orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=8,dimension=3)
orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=8,dimension=4)
orderlist1_AB5<-createorderlist(quni1=quasiuni_sorted5,size=samplesize,interval=8,dimension=5)

batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)



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
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted = TRUE)
    midhinge1<-midhinge(x=sortedx,sorted = TRUE)
    SWA8<-SWA(x=sortedx,interval=8,batch="auto",sorted = TRUE)
    SWAHlmean1<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,orderlist1_sorted5=orderlist1_AB5,interval=8,batch="auto")
    THL1<-THL(x=sortedx,orderlist1_sorted2=orderlist1_AB2,interval=4,batch="auto",boot=TRUE)
    rqm1<-rqm(x=sortedx,sorted = TRUE)
    MoM2<-median_of_means(x=sortedx,korder=2)
    MoM3<-median_of_means(x=sortedx,korder=3)
    MoM4<-median_of_means(x=sortedx,korder=4)
    MoM5<-median_of_means(x=sortedx,korder=5)
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,midhinge1,SWA8,SWAHlmean1,THL1,rqm1,MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5)-momentsx[1])/(sqrt(momentsx[2]))
    )
    
    all1<-t(c(kurtx,Huberx,SMWM9,midhinge1,SWA8,SWAHlmean1,THL1,rqm1,MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,targetall,momentsx,allrawmoBias))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("gnorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,61:111])^2))/sqrt(targetvar)
  
  AB1_mean<-abs(colMeans((SEbataches[,61:111])))/sqrt(targetvar)
  
  SEbatachesmean <- colMeans(SEbataches)
  
  ratiomean1<-c(SEbatachesmean[c(1:60)])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(1:60)]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(1:60)])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],RMSE1_mean=RMSE1_mean,AB1_mean=AB1_mean,SEbatachesmean=SEbatachesmean,meansd_unscaled1=meansd_unscaled1,mean_SE1=mean_SE1,meansd1=meansd1,mean_SSE1=mean_SSE1)

  allErrors
}

write.csv(simulatedbatch_ABSE,paste("gnorm_ABSSE.csv", sep = ","), row.names = FALSE)

