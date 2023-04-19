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

numCores <- detectCores()-4 # Detect the number of available cores
cl <- makeCluster(numCores) # Create a cluster with the number of cores
registerDoParallel(cl) # Register the parallel backend


#bootsize for bootstrap approximation of the distributions of the kernal of U-statistics.
n <- 13824*2*3
(n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (approximating asymptotic)
largesize<-13824*2

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)

quasiuni<-quasiunisobol

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))

# Forever...

#load asymptotic d for two parameter distributions

#set the stop criterion
criterionset=1e-10

kurtWeibull<- read.csv(("kurtWeibull_28260.csv"))

allkurtWeibull<-unlist(kurtWeibull)

samplesize=576*9
batchsizebase=1000

orderlist1_AB20<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=8,dimension=2)
orderlist1_AB20<-orderlist1_AB20[1:largesize,]
orderlist1_AB30<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=8,dimension=3)
orderlist1_AB30<-orderlist1_AB30[1:largesize,]
orderlist1_AB40<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=8,dimension=4)
orderlist1_AB40<-orderlist1_AB40[1:largesize,]

batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

#input the d value table previously generated
d_values<- read.csv(("d_values.csv"))
I_values<-read.csv(("I_values.csv"))
Ismoments_values<-read.csv(("Ismoments_values.csv"))
#Then, start the Monte Simulation

simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  set.seed(1)
  a=allkurtWeibull[batchnumber]
  
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-c(kurtx=targetfm/(targetvar^(4/2)))
  skewx<-c(skewx=targettm/(targetvar^(3/2)))
  
  SEbataches<- read.csv(paste("finite_Weibull_Icalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  SEbataches2<- read.csv(paste("finite_Weibull_Ismomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbataches3<-c()
  for (batch1 in c(1:batchsize)){
    iall11<-SEbataches[batch1,]
    iall12<-SEbataches2[batch1,]
    
    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    
    rqmomentselect1<-rqmoments3(x=sortedx,iall1=iall11,ismoments1=iall12,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_Ismoments=Ismoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE)
    
    standardizedmomentsx<-standardizedmoments(x=sortedx)
    
    sortedx<-c()
    
    all1<-(c(rqmomentselect1,targetall,standardizedmomentsx))
    
    SEbataches3<-rbind(SEbataches3,all1)
  }
  
  write.csv(SEbataches3,paste("finite_Weibull_Imomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  SEbatachesmean <-apply(SEbataches3, 2, calculate_column_mean)
  
  rqmean<-apply(((SEbataches3[1:batchsize,c(3:74,195:266,387:458,579:650)])), 2, calculate_column_sd)
  
  rqvar<-apply((SEbataches3[1:batchsize,c(75:126,267:318,459:510,651:702)]), 2, calculate_column_sd)
  
  rqtm<-apply((SEbataches3[1:batchsize,c(127:166,319:358,511:550,703:742)]), 2, calculate_column_sd)
  
  rqfm<-apply((SEbataches3[1:batchsize,c(167:194,359:386,551:578,743:70)]), 2, calculate_column_sd)
  
  rankmean1<-rank(rqmean)
  rankvar1<-rank(rqvar)
  ranktm1<-rank(rqtm)
  rankfm1<-rank(rqfm)
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,rankmean1,rankvar1,ranktm1,rankfm1,SEbatachesmean,SErqmean=rqmean,SErqvar=rqvar,SErqtm=rqtm,SErqfm=rqfm)
}


write.csv(simulatedbatch_bias_Monte,paste("finite_Weibull_Imomentscalibration_raw",samplesize,".csv", sep = ","), row.names = FALSE)

Optimum_SE<-simulatedbatch_bias_Monte[,1:1418]

write.csv(Optimum_SE,paste("finite_Imoments_Weibull.csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
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
  
  SEbataches<- read.csv(paste("finite_Weibull_Imomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  
  
  rqmean_se<-apply(((SEbataches[1:batchsize,c(3:74,195:266,387:458,579:650)])), 2, se_sd)
  
  rqvar_se<-apply((SEbataches[1:batchsize,c(75:126,267:318,459:510,651:702)]), 2, se_sd)
  
  rqtm_se<-apply((SEbataches[1:batchsize,c(127:166,319:358,511:550,703:742)]), 2, se_sd)
  
  rqfm_se<-apply((SEbataches[1:batchsize,c(167:194,359:386,551:578,743:70)]), 2, se_sd)
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,se_mean_all1,rqmean_se,rqvar_se,rqtm_se,rqfm_se)
  
  allresultsSE
}

write.csv(simulatedbatch_bias_Monte_SE,paste("finite_Weibull_Imomentscalibration_raw_error",samplesize,".csv", sep = ","), row.names = FALSE)

finite_I_Weibull<- read.csv(("finite_Imoments_Weibull.csv"))


finite_I_Weibull<-data.frame(finite_I_Weibull)

asymptotic_I_Weibull<- read.csv(paste("asymptotic_Imoments.csv", sep = ","))
colnames(asymptotic_I_Weibull)<-colnames(finite_I_Weibull)

all1<-rbind(asymptotic_I_Weibull,finite_I_Weibull)

write.csv(all1,paste("Imoments_values.csv", sep = ","), row.names = FALSE)

stopCluster(cl)
registerDoSEQ()

