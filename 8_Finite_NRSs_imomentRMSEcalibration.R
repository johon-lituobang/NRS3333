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

kurtlognorm<- read.csv(("kurtlognorm_31260.csv"))
allkurtlognorm<-unlist(kurtlognorm)


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

simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:length(allkurtlognorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  set.seed(1)
  a=allkurtlognorm[batchnumber]
  
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  
  RMSEbataches<- read.csv(paste("finite_lognorm_Icalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  RMSEbataches2<- read.csv(paste("finite_lognorm_Ismomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  RMSEbataches3<-c()
  for (batch1 in c(1:batchsize)){
    iall11<-RMSEbataches[batch1,]
    iall12<-RMSEbataches2[batch1,]
    
    x<-c(dslnorm(uni=unibatch[,batch1], meanlog =0, sdlog  = a/1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    
    rqmomentselect1<-rqmoments3(x=sortedx,iall1=iall11,ismoments1=iall12,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_Ismoments=Ismoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE)
    
    standardizedmomentsx<-standardizedmoments(x=sortedx)
    
    sortedx<-c()
    
    all1<-(c(rqmomentselect1,targetall,standardizedmomentsx))
    
    RMSEbataches3<-rbind(RMSEbataches3,all1)
  }
  
  write.csv(RMSEbataches3,paste("finite_lognorm_Imomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSEbatachesmean <-apply(RMSEbataches3, 2, calculate_column_mean)
  
  rqmean<-sqrt(colMeans((RMSEbataches3[1:batchsize,c(3:74,195:266,387:458,579:650,771:842,963:1034,1155:1226,1347:1418,1539:1610,1731:1802)]-targetm)^2))
  
  rqvar<-sqrt(colMeans((RMSEbataches3[1:batchsize,c(75:126,267:318,459:510,651:702,843:894,1035:1086,1227:1278,1419:1470,1611:1662,1803:1854)]-targetvar)^2))
  
  rqtm<-sqrt(colMeans((RMSEbataches3[1:batchsize,c(127:166,319:358,511:550,703:742,895:934,1087:1126,1279:1318,1471:1510,1663:1702,1855:1894)]-targettm)^2))
  
  rqfm<-sqrt(colMeans((RMSEbataches3[1:batchsize,c(167:194,359:386,551:578,743:770,935:962,1127:1154,1319:1346,1511:1538,1703:1730,1895:1922)]-targetfm)^2))
  
  rankmean1<-rank(rqmean)
  rankvar1<-rank(rqvar)
  ranktm1<-rank(rqtm)
  rankfm1<-rank(rqfm)
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,rankmean1,rankvar1,ranktm1,rankfm1,RMSEbatachesmean,RMSErqmean=rqmean,RMSErqvar=rqvar,RMSErqtm=rqtm,RMSErqfm=rqfm)
}


write.csv(simulatedbatch_bias_Monte,paste("finite_lognorm_Imomentscalibration_raw",samplesize,".csv", sep = ","), row.names = FALSE)

Optimum_RMSE<-simulatedbatch_bias_Monte[,1:1924]

write.csv(Optimum_RMSE,paste("finite_Imoments_lognorm.csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:length(allkurtlognorm))), .combine = 'rbind') %dopar% {
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
  
  
  SEbataches<- read.csv(paste("finite_lognorm_Imomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  
  
  rqmean_se<-apply(((SEbataches[1:batchsize,c(3:74,195:266,387:458,579:650,771:842,963:1034,1155:1226,1347:1418,1539:1610,1731:1802)])), 2, se_sd)
  
  rqvar_se<-apply((SEbataches[1:batchsize,c(75:126,267:318,459:510,651:702,843:894,1035:1086,1227:1278,1419:1470,1611:1662,1803:1854)]), 2, se_sd)
  
  rqtm_se<-apply((SEbataches[1:batchsize,c(127:166,319:358,511:550,703:742,895:934,1087:1126,1279:1318,1471:1510,1663:1702,1855:1894)]), 2, se_sd)
  
  rqfm_se<-apply((SEbataches[1:batchsize,c(167:194,359:386,551:578,743:770,935:962,1127:1154,1319:1346,1511:1538,1703:1730,1895:1922)]), 2, se_sd)
  
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,se_mean_all1,rqmean_se,rqvar_se,rqtm_se,rqfm_se)
  
  allresultsSE
}

write.csv(simulatedbatch_bias_Monte_SE,paste("finite_lognorm_Imomentscalibration_raw_error",samplesize,".csv", sep = ","), row.names = FALSE)

finite_I_lognorm<- read.csv(("finite_Imoments_lognorm.csv"))


finite_I_lognorm<-data.frame(finite_I_lognorm)

asymptotic_I_lognorm<- read.csv(paste("asymptotic_Imoments.csv", sep = ","))
colnames(asymptotic_I_lognorm)<-colnames(finite_I_lognorm)

all1<-rbind(asymptotic_I_lognorm,finite_I_lognorm)

write.csv(all1,paste("Imoments_values.csv", sep = ","), row.names = FALSE)

stopCluster(cl)
registerDoSEQ()

