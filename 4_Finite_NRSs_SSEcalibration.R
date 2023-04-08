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
n <- 2048*9*3
(n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (approximating asymptotic)
largesize<-2048*9

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

samplesize=2048*2
batchsizebase=1000

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
batchsize=1000

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

#input the d value table previously generated
d_values<- read.csv(("d_SWA.csv"))
#Then, start the Monte Simulation

simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  setSeed(batchnumber)
  set.seed(batchnumber)
  a=allkurtWeibull[batchnumber]
  
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-c(kurtx=targetfm/(targetvar^(4/2)))
  skewx<-c(skewx=targettm/(targetvar^(3/2)))
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    
    rqmoments1<-rqsmoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset)
    
    standardizedmomentsx<-standardizedmoments(x=sortedx)
    
    sortedx<-c()
    
    all1<-t(c(rqmoments1,targetall,standardizedmomentsx))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("finite_Weibull_Icalibration_raw_SWA",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  SEbatachesmean <-apply(SEbataches, 2, calculate_column_mean)
  
  rqkurt<-apply(((SEbataches[1:batchsize,c(1:238)])), 2, calculate_column_sd)
  
  rqskew<-apply((SEbataches[1:batchsize,c(239:578)]), 2, calculate_column_sd)
  
  rankkurtall1<-rank(rqkurt[c(1:238)])
  rankskewall1<-rank(rqskew[c(1:340)])
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,rankkurtall1,rankskewall1,SEbatachesmean,SErqkurt=rqkurt,SErqskew=rqskew)
}

write.csv(simulatedbatch_bias_Monte,paste("finite_Weibull_Icalibration_raw_SWA",samplesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte<- read.csv(paste("finite_Weibull_Icalibration_raw_SWA",samplesize,".csv", sep = ","))

Optimum_SE<-simulatedbatch_bias_Monte[,1:582]

write.csv(Optimum_SE,paste("finite_I_Weibull_SWA.csv", sep = ","), row.names = FALSE)

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
  
  SEbataches<- read.csv(paste("finite_Weibull_Icalibration_raw_SWA",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  
  rqkurt_se<-apply(((SEbataches[1:batchsize,c(1:238)])), 2, se_sd)
  
  rqskew_se<-apply((SEbataches[1:batchsize,c(239:578)]), 2, se_sd)
  
  allresultsSE<-c(samplesize=samplesize,type=1,kurtx,skewx,se_mean_all1,rqkurt_se,rqskew_se)
  
  allresultsSE
}

write.csv(simulatedbatch_bias_Monte_SE,paste("finite_Weibull_Icalibration_raw_SWA_error",samplesize,".csv", sep = ","), row.names = FALSE)


finite_I_Weibull<- read.csv(("finite_I_Weibull_SWA.csv"))


finite_I_Weibull<-data.frame(finite_I_Weibull)

asymptotic_I_Weibull<- read.csv(paste("asymptotic_I_SWA.csv", sep = ","))
colnames(finite_I_Weibull)<-colnames(asymptotic_I_Weibull)
all1<-rbind(asymptotic_I_Weibull,finite_I_Weibull)

write.csv(all1,paste("I_SWA.csv", sep = ","), row.names = FALSE)



stopCluster(cl)
registerDoSEQ()

