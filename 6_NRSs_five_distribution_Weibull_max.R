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
if (!require("NRSReview")) install.packages("NRSReview_1.0.zip", repos = NULL)
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

simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:100), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=0.01*batchnumber
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
  imoments1<-tryCatch({
    imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,171)
  })
  
  momentsx<-unbiasedmoments(x=sortedx)

  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-momentsx[1])/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-momentsx[2])/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-momentsx[3])/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-momentsx[4])/momentssd[4])
  
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
}
nametwo<- read.csv(paste("asymptotic_max_raw_Process_label.csv", sep = ","))
colnames(simulatedbatch_asymptoticbias)<-colnames(nametwo)
write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_max_Weibull_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:100,1],simulatedbatch_asymptoticbias[1:100,193:375]),paste("asymptotic_max_Weibull",largesize,".csv", sep = ","), row.names = FALSE)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:100), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=0.01*batchnumber
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
  imoments1<-tryCatch({
    imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,171)
  })
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-momentsx[1])/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-momentsx[2])/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-momentsx[3])/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-momentsx[4])/momentssd[4])
  
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
}
nametwo<- read.csv(paste("asymptotic_max_raw_Process_label.csv", sep = ","))
colnames(simulatedbatch_asymptoticbias)<-colnames(nametwo)
write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_max_gamma_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:100,1],simulatedbatch_asymptoticbias[1:100,193:375]),paste("asymptotic_max_gamma",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:100), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=4+0.1*batchnumber
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
  imoments1<-tryCatch({
    imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,171)
  })
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-momentsx[1])/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-momentsx[2])/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-momentsx[3])/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-momentsx[4])/momentssd[4])
  
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
}
nametwo<- read.csv(paste("asymptotic_max_raw_Process_label.csv", sep = ","))
colnames(simulatedbatch_asymptoticbias)<-colnames(nametwo)
write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_max_Pareto_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:100,1],simulatedbatch_asymptoticbias[1:100,193:375]),paste("asymptotic_max_Pareto",largesize,".csv", sep = ","), row.names = FALSE)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:100), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=batchnumber
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
  imoments1<-tryCatch({
    imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,171)
  })
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-momentsx[1])/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-momentsx[2])/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-momentsx[3])/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-momentsx[4])/momentssd[4])
  
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
}
nametwo<- read.csv(paste("asymptotic_max_raw_Process_label.csv", sep = ","))
colnames(simulatedbatch_asymptoticbias)<-colnames(nametwo)
write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_max_lognorm_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:100,1],simulatedbatch_asymptoticbias[1:100,193:375]),paste("asymptotic_max_lognorm",largesize,".csv", sep = ","), row.names = FALSE)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:100), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  a=0.01*batchnumber
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
  imoments1<-tryCatch({
    imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,171)
  })
  momentsx<-unbiasedmoments(x=sortedx)
  
  sortedx<-c()
  
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[153],imoments1[154],imoments1[155])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,HL1,imoments1[seq(from=1, to=96, by=4)],imoments1[97:110],imoments1[seq(from=156, to=171, by=4)],momentsx[1])-momentsx[1])/momentssd[1],
    secondbias=abs(c(imoments1[seq(from=2, to=96, by=4)],imoments1[111:124],imoments1[seq(from=157, to=171, by=4)],momentsx[2])-momentsx[2])/momentssd[2],
    thirdbias=abs(c(imoments1[seq(from=3, to=96, by=4)],imoments1[125:138],imoments1[seq(from=158, to=171, by=4)],momentsx[3])-momentsx[3])/momentssd[3],
    fourbias=abs(c(imoments1[seq(from=4, to=96, by=4)],imoments1[139:152],imoments1[seq(from=159, to=171, by=4)],momentsx[4])-momentsx[4])/momentssd[4])
  
  all1<-t(c(kurtx,skewx,Huberx,SMWM9,HL1,imoments1,targetall,momentsx,allrawmoBias,momentssd))
}
nametwo<- read.csv(paste("asymptotic_max_raw_Process_label.csv", sep = ","))
colnames(simulatedbatch_asymptoticbias)<-colnames(nametwo)
write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_max_gnorm_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:100,1],simulatedbatch_asymptoticbias[1:100,193:375]),paste("asymptotic_max_gnorm",largesize,".csv", sep = ","), row.names = FALSE)

