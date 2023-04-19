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
n <- 331776*3*8
(n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (approximating asymptotic)
largesize<-331776*8

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

samplesize=331776*8
batchsizebase=5


orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2,size=largesize,interval=8,dimension=2)
orderlist1_AB2<-orderlist1_AB2[1:largesize,]
orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3,size=largesize,interval=8,dimension=3)
orderlist1_AB3<-orderlist1_AB3[1:largesize,]
orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4,size=largesize,interval=8,dimension=4)
orderlist1_AB4<-orderlist1_AB4[1:largesize,]
batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

#input the d value table previously generated

#Then, start the Monte Simulation

#input the d value table previously generated
d_values<- read.csv(("d_values.csv"))
I_values<-read.csv(("I_values.csv"))
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
  SEbataches<- read.csv(paste("asymptotic_Weibull_Icalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbataches2<-c()
  for (batch1 in c(1:batchsize)){
    iall11<-SEbataches[batch1,]
    
    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
      
    rqmomentselect1<-rqmoments2(x=sortedx,iall1=iall11,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted20=orderlist1_AB2,orderlist1_sorted30=orderlist1_AB3,orderlist1_sorted40=orderlist1_AB4,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE)
    
    standardizedmomentsx<-standardizedmoments(x=sortedx)
    
    sortedx<-c()
    
    all1<-(c(rqmomentselect1,targetall,standardizedmomentsx))
    
    SEbataches2<-rbind(SEbataches2,all1)
  }
  
  write.csv(SEbataches2,paste("asymptotic_Weibull_Imomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  SEbatachesmean <-apply(SEbataches2, 2, calculate_column_mean)
  
  rqmean<-apply(((SEbataches2[1:batchsize,1:72])), 2, calculate_column_sd)
  
  rqvar<-apply((SEbataches2[1:batchsize,c(73:124)]), 2, calculate_column_sd)
  
  rqtm<-apply((SEbataches2[1:batchsize,c(125:164)]), 2, calculate_column_sd)
  
  rqfm<-apply((SEbataches2[1:batchsize,c(165:192)]), 2, calculate_column_sd)
  
  rankmean1<-rank(rqmean)
  rankvar1<-rank(rqvar)
  ranktm1<-rank(rqtm)
  rankfm1<-rank(rqfm)
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,rankmean1,rankvar1,ranktm1,rankfm1,SEbatachesmean,SErqmean=rqmean,SErqvar=rqvar,SErqtm=rqtm,SErqfm=rqfm)
}

write.csv(simulatedbatch_bias_Monte,paste("asymptotic_Weibull_Imomentscalibration_raw",largesize,".csv", sep = ","), row.names = FALSE)

Optimum_SE<-simulatedbatch_bias_Monte[,1:196]

write.csv(Optimum_SE,paste("asymptotic_Imoments_Weibull.csv", sep = ","), row.names = FALSE)

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

  SEbataches<- read.csv(paste("asymptotic_Weibull_Imomentscalibration_raw",samplesize,round(kurtx,digits = 1),".csv", sep = ","))

  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  
  rqmean_se<-apply(((SEbataches[1:batchsize,1:72])), 2, se_sd)
  
  rqvar_se<-apply((SEbataches[1:batchsize,c(73:124)]), 2, se_sd)
  
  rqtm_se<-apply((SEbataches[1:batchsize,c(125:164)]), 2, se_sd)
  
  rqfm_se<-apply((SEbataches[1:batchsize,c(165:192)]), 2, se_sd)
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,se_mean_all1,rqmean_se,rqvar_se,rqtm_se,rqfm_se)

  allresultsSE
}

write.csv(simulatedbatch_bias_Monte_SE,paste("asymptotic_Weibull_Imomentscalibration_raw_error",largesize,".csv", sep = ","), row.names = FALSE)

asymptotic_I_Weibull<- read.csv(("asymptotic_Imoments_Weibull.csv"))


write.csv(asymptotic_I_Weibull,paste("asymptotic_Imoments.csv", sep = ","), row.names = FALSE)

stopCluster(cl)
registerDoSEQ()

