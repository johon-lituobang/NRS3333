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



asymptotic_n <- 331776*3*8
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (asymptotic bias)
largesize<-331776*8

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_asymptotic<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

quasiuni_sorted2_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4_asymptotic <- na.omit(rowSort(quasiuni_asymptotic, descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...
quasiuni_asymptotic<-Sort(sobol(n=largesize, dim = 1, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1))


orderlist1_AB2_asymptotic<-createorderlist(quni1=quasiuni_sorted2_asymptotic,size=largesize,interval=16,dimension=2)
orderlist1_AB2_asymptotic<-orderlist1_AB2_asymptotic[1:largesize,]
orderlist1_AB3_asymptotic<-createorderlist(quni1=quasiuni_sorted3_asymptotic,size=largesize,interval=16,dimension=3)
orderlist1_AB3_asymptotic<-orderlist1_AB3_asymptotic[1:largesize,]
orderlist1_AB4_asymptotic<-createorderlist(quni1=quasiuni_sorted4_asymptotic,size=largesize,interval=16,dimension=4)
orderlist1_AB4_asymptotic<-orderlist1_AB4_asymptotic[1:largesize,]

quasiuni_sorted2_asymptotic<-c()
quasiuni_sorted3_asymptotic<-c()
quasiuni_sorted4_asymptotic<-c()


asymptotic_n <- 331776*3*8
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 6
#large sample size (asymptotic bias)
largesize<-331776*8

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_M<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

orderlist1_hllarge_asymptotic<-createorderlist(quni1=quasiuni_M[,1:6],size=largesize,interval=8,dimension=6)
orderlist1_hllarge_asymptotic<-orderlist1_hllarge_asymptotic[1:largesize,]

morder=6
unibatchran_M<-matrix(randtoolbox::SFMT(largesize*3*morder),ncol=morder)

orderlist1_hllarge_asymptotic_rand<-createorderlist(quni1=unibatchran_M[,1:6],size=largesize,interval=8,dimension=6)
orderlist1_hllarge_asymptotic_rand<-orderlist1_hllarge_asymptotic_rand[1:largesize,]


#load asymptotic d for two parameter distributions
d_values<- read.csv(("d_values.csv"))
I_values<- read.csv(("I_values.csv"))
Ismoments_values<- read.csv(("Ismoments_values.csv"))
Imoments_values<- read.csv(("Imoments_values.csv"))

#set the stop criterion
criterionset=1e-10

kurtWeibull<- read.csv(("kurtWeibull_28100000.csv"))
allkurtWeibull<-unlist(kurtWeibull)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtWeibull)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
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
  
  Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  
  imoments1<-tryCatch({
    unlist(imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Ismoments=Ismoments_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_hlsmall=orderlist1_hllarge_asymptotic,orderlist1_hllarge=orderlist1_hllarge_asymptotic,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE))
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,3903)
  })
  momentsx<-standardizedmoments(x=sortedx,releaseall = TRUE)
  
  #D Olive, Robust estimators for transformed location scale families. Unpubl. manuscript 1025 available from (www. math. siu. edu/olive/preprints. htm) (2006).
  medianMAD1<-Weibull_median_MAD_estimator(sortedx)
  
  #NB Marks, Estimation of weibull parameters from common percentiles. J. applied Stat. 32, 17–24 (2005). 
  QE1<-Weibull_quantile_estimator(x=sortedx,sorted=TRUE)
  
  alpha1<-QE1[1]-0.3
  alpha2<-QE1[1]+0.3
  
  #X He, WK Fung, Method of medians for lifetime data with weibull models. Stat. medicine 18, 1993–2009 (1999)
  RMLE1<-Weibull_RMLE(sortedx,alpha1=alpha1,alpha2=alpha2)
  
  #all parameter setting are from
  #K Boudt, D Caliskan, C Croux, Robust explicit estimators of weibull parameters. Metrika 73, 187–209 (2011).
  
  moments_medianMAD1<-Weibull_moments(alpha=medianMAD1[1],lambda=medianMAD1[2])
  
  moments_QE1<-Weibull_moments(alpha=QE1[1],lambda=QE1[2])
  
  moments_RMLE1<-Weibull_moments(alpha=RMLE1[1],lambda=RMLE1[2])
  
  MoM2<-median_of_means(sortedx,korder=2)
  MoM3<-median_of_means(sortedx,korder=3)
  MoM4<-median_of_means(sortedx,korder=4)
  MoM5<-median_of_means(sortedx,korder=5)
  MoM6<-median_of_means(sortedx,korder=6)
  
  MoRMall<-mHLM_all(x=sortedx,max_dim=6,orderlists=orderlist1_hllarge_asymptotic_rand,boot=TRUE,quasi=FALSE,largesize=largesize)
  
  MoRMall<-unlist(MoRMall)
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[3899],imoments1[3900],imoments1[3901])

  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)-imoments1[3699])/momentssd[1],
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-imoments1[3773])/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-imoments1[3827])/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-imoments1[3869])/momentssd[4],
    skewbias=abs(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])-imoments1[3903]),
    kurtbias=abs(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902])-imoments1[3902])
  )
  
  allrawmo1<-c(first=(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)),
               second=(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])),
               third=(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])),
               fourth=(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])),
               skew=(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])),
               kurt=(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902]))
  )
  
  
  medianmoments<-c(imoments1[3706],imoments1[3780],imoments1[3834],imoments1[3875])
  standardizedm<-c(imoments1[3706]/momentssd[1],imoments1[3780]/momentssd[2],imoments1[3834]/momentssd[3],imoments1[3875]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_Weibull_raw_Process_max",largesize,".csv", sep = ","), row.names = FALSE)


kurtgamma<- read.csv(("kurtgamma_31100000.csv"))
allkurtgamma<-unlist(kurtgamma)

simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtgamma)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
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
  
  Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  
  imoments1<-tryCatch({
    unlist(imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Ismoments=Ismoments_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_hlsmall=orderlist1_hllarge_asymptotic,orderlist1_hllarge=orderlist1_hllarge_asymptotic,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE))
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,3903)
  })
  momentsx<-standardizedmoments(x=sortedx,releaseall = TRUE)
  
  #D Olive, Robust estimators for transformed location scale families. Unpubl. manuscript 1025 available from (www. math. siu. edu/olive/preprints. htm) (2006).
  medianMAD1<-Weibull_median_MAD_estimator(sortedx)
  
  #NB Marks, Estimation of weibull parameters from common percentiles. J. applied Stat. 32, 17–24 (2005). 
  QE1<-Weibull_quantile_estimator(x=sortedx,sorted=TRUE)
  
  alpha1<-QE1[1]-0.3
  alpha2<-QE1[1]+0.3
  
  #X He, WK Fung, Method of medians for lifetime data with weibull models. Stat. medicine 18, 1993–2009 (1999)
  RMLE1<-Weibull_RMLE(sortedx,alpha1=alpha1,alpha2=alpha2)
  
  #all parameter setting are from
  #K Boudt, D Caliskan, C Croux, Robust explicit estimators of weibull parameters. Metrika 73, 187–209 (2011).
  
  moments_medianMAD1<-Weibull_moments(alpha=medianMAD1[1],lambda=medianMAD1[2])
  
  moments_QE1<-Weibull_moments(alpha=QE1[1],lambda=QE1[2])
  
  moments_RMLE1<-Weibull_moments(alpha=RMLE1[1],lambda=RMLE1[2])
  
  MoM2<-median_of_means(sortedx,korder=2)
  MoM3<-median_of_means(sortedx,korder=3)
  MoM4<-median_of_means(sortedx,korder=4)
  MoM5<-median_of_means(sortedx,korder=5)
  MoM6<-median_of_means(sortedx,korder=6)
  
  MoRMall<-mHLM_all(x=sortedx,max_dim=6,orderlists=orderlist1_hllarge_asymptotic_rand,boot=TRUE,quasi=FALSE,largesize=largesize)
  
  MoRMall<-unlist(MoRMall)
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[3899],imoments1[3900],imoments1[3901])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)-imoments1[3699])/momentssd[1],
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-imoments1[3773])/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-imoments1[3827])/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-imoments1[3869])/momentssd[4],
    skewbias=abs(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])-imoments1[3903]),
    kurtbias=abs(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902])-imoments1[3902])
  )
  
  allrawmo1<-c(first=(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)),
               second=(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])),
               third=(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])),
               fourth=(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])),
               skew=(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])),
               kurt=(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902]))
  )
  
  
  medianmoments<-c(imoments1[3706],imoments1[3780],imoments1[3834],imoments1[3875])
  standardizedm<-c(imoments1[3706]/momentssd[1],imoments1[3780]/momentssd[2],imoments1[3834]/momentssd[3],imoments1[3875]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_gamma_raw_Process_max",largesize,".csv", sep = ","), row.names = FALSE)


kurtPareto<- read.csv(("kurtPareto_91100000.csv"))
allkurtPareto<-unlist(kurtPareto)

simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtPareto)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
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
  
  
  Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  
  imoments1<-tryCatch({
    unlist(imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Ismoments=Ismoments_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_hlsmall=orderlist1_hllarge_asymptotic,orderlist1_hllarge=orderlist1_hllarge_asymptotic,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE))
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,3903)
  })
  momentsx<-standardizedmoments(x=sortedx,releaseall = TRUE)
  
  #D Olive, Robust estimators for transformed location scale families. Unpubl. manuscript 1025 available from (www. math. siu. edu/olive/preprints. htm) (2006).
  medianMAD1<-Weibull_median_MAD_estimator(sortedx)
  
  #NB Marks, Estimation of weibull parameters from common percentiles. J. applied Stat. 32, 17–24 (2005). 
  QE1<-Weibull_quantile_estimator(x=sortedx,sorted=TRUE)
  
  alpha1<-QE1[1]-0.3
  alpha2<-QE1[1]+0.3
  
  #X He, WK Fung, Method of medians for lifetime data with weibull models. Stat. medicine 18, 1993–2009 (1999)
  RMLE1<-Weibull_RMLE(sortedx,alpha1=alpha1,alpha2=alpha2)
  
  #all parameter setting are from
  #K Boudt, D Caliskan, C Croux, Robust explicit estimators of weibull parameters. Metrika 73, 187–209 (2011).
  
  moments_medianMAD1<-Weibull_moments(alpha=medianMAD1[1],lambda=medianMAD1[2])
  
  moments_QE1<-Weibull_moments(alpha=QE1[1],lambda=QE1[2])
  
  moments_RMLE1<-Weibull_moments(alpha=RMLE1[1],lambda=RMLE1[2])
  
  MoM2<-median_of_means(sortedx,korder=2)
  MoM3<-median_of_means(sortedx,korder=3)
  MoM4<-median_of_means(sortedx,korder=4)
  MoM5<-median_of_means(sortedx,korder=5)
  MoM6<-median_of_means(sortedx,korder=6)
  
  MoRMall<-mHLM_all(x=sortedx,max_dim=6,orderlists=orderlist1_hllarge_asymptotic_rand,boot=TRUE,quasi=FALSE,largesize=largesize)
  
  MoRMall<-unlist(MoRMall)
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[3899],imoments1[3900],imoments1[3901])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)-imoments1[3699])/momentssd[1],
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-imoments1[3773])/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-imoments1[3827])/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-imoments1[3869])/momentssd[4],
    skewbias=abs(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])-imoments1[3903]),
    kurtbias=abs(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902])-imoments1[3902])
  )
  
  allrawmo1<-c(first=(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)),
               second=(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])),
               third=(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])),
               fourth=(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])),
               skew=(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])),
               kurt=(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902]))
  )
  
  
  medianmoments<-c(imoments1[3706],imoments1[3780],imoments1[3834],imoments1[3875])
  standardizedm<-c(imoments1[3706]/momentssd[1],imoments1[3780]/momentssd[2],imoments1[3834]/momentssd[3],imoments1[3875]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_Pareto_raw_Process_max",largesize,".csv", sep = ","), row.names = FALSE)


kurtlognorm<- read.csv(("kurtlognorm_31100000.csv"))
allkurtlognorm<-unlist(kurtlognorm)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtlognorm)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
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
  
  
  Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  
  imoments1<-tryCatch({
    unlist(imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Ismoments=Ismoments_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_hlsmall=orderlist1_hllarge_asymptotic,orderlist1_hllarge=orderlist1_hllarge_asymptotic,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE))
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,3903)
  })
  momentsx<-standardizedmoments(x=sortedx,releaseall = TRUE)
  
  #D Olive, Robust estimators for transformed location scale families. Unpubl. manuscript 1025 available from (www. math. siu. edu/olive/preprints. htm) (2006).
  medianMAD1<-Weibull_median_MAD_estimator(sortedx)
  
  #NB Marks, Estimation of weibull parameters from common percentiles. J. applied Stat. 32, 17–24 (2005). 
  QE1<-Weibull_quantile_estimator(x=sortedx,sorted=TRUE)
  
  alpha1<-QE1[1]-0.3
  alpha2<-QE1[1]+0.3
  
  #X He, WK Fung, Method of medians for lifetime data with weibull models. Stat. medicine 18, 1993–2009 (1999)
  RMLE1<-Weibull_RMLE(sortedx,alpha1=alpha1,alpha2=alpha2)
  
  #all parameter setting are from
  #K Boudt, D Caliskan, C Croux, Robust explicit estimators of weibull parameters. Metrika 73, 187–209 (2011).
  
  moments_medianMAD1<-Weibull_moments(alpha=medianMAD1[1],lambda=medianMAD1[2])
  
  moments_QE1<-Weibull_moments(alpha=QE1[1],lambda=QE1[2])
  
  moments_RMLE1<-Weibull_moments(alpha=RMLE1[1],lambda=RMLE1[2])
  
  MoM2<-median_of_means(sortedx,korder=2)
  MoM3<-median_of_means(sortedx,korder=3)
  MoM4<-median_of_means(sortedx,korder=4)
  MoM5<-median_of_means(sortedx,korder=5)
  MoM6<-median_of_means(sortedx,korder=6)
  
  MoRMall<-mHLM_all(x=sortedx,max_dim=6,orderlists=orderlist1_hllarge_asymptotic_rand,boot=TRUE,quasi=FALSE,largesize=largesize)
  
  MoRMall<-unlist(MoRMall)
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[3899],imoments1[3900],imoments1[3901])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)-imoments1[3699])/momentssd[1],
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-imoments1[3773])/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-imoments1[3827])/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-imoments1[3869])/momentssd[4],
    skewbias=abs(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])-imoments1[3903]),
    kurtbias=abs(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902])-imoments1[3902])
  )
  
  allrawmo1<-c(first=(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)),
               second=(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])),
               third=(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])),
               fourth=(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])),
               skew=(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])),
               kurt=(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902]))
  )
  
  
  medianmoments<-c(imoments1[3706],imoments1[3780],imoments1[3834],imoments1[3875])
  standardizedm<-c(imoments1[3706]/momentssd[1],imoments1[3780]/momentssd[2],imoments1[3834]/momentssd[3],imoments1[3875]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_lognorm_raw_Process_max",largesize,".csv", sep = ","), row.names = FALSE)


kurtgnorm<- read.csv(("kurtgnorm_21100000.csv"))
allkurtgnorm<-unlist(kurtgnorm)


simulatedbatch_asymptoticbias<-foreach(batchnumber = (1:length(allkurtgnorm)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
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
  
  
  Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
  
  imoments1<-tryCatch({
    unlist(imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Ismoments=Ismoments_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_hlsmall=orderlist1_hllarge_asymptotic,orderlist1_hllarge=orderlist1_hllarge_asymptotic,percentage=1/24,batch="auto",stepsize=1000,criterion=1e-10,boot=TRUE))
  }, error = function(e) {
    cat("Error: ", conditionMessage(e), "\n")
    rep(NA,3903)
  })
  momentsx<-standardizedmoments(x=sortedx,releaseall = TRUE)
  
  #D Olive, Robust estimators for transformed location scale families. Unpubl. manuscript 1025 available from (www. math. siu. edu/olive/preprints. htm) (2006).
  medianMAD1<-Weibull_median_MAD_estimator(sortedx)
  
  #NB Marks, Estimation of weibull parameters from common percentiles. J. applied Stat. 32, 17–24 (2005). 
  QE1<-Weibull_quantile_estimator(x=sortedx,sorted=TRUE)
  
  alpha1<-QE1[1]-0.3
  alpha2<-QE1[1]+0.3
  
  #X He, WK Fung, Method of medians for lifetime data with weibull models. Stat. medicine 18, 1993–2009 (1999)
  RMLE1<-Weibull_RMLE(sortedx,alpha1=alpha1,alpha2=alpha2)
  
  #all parameter setting are from
  #K Boudt, D Caliskan, C Croux, Robust explicit estimators of weibull parameters. Metrika 73, 187–209 (2011).
  
  moments_medianMAD1<-Weibull_moments(alpha=medianMAD1[1],lambda=medianMAD1[2])
  
  moments_QE1<-Weibull_moments(alpha=QE1[1],lambda=QE1[2])
  
  moments_RMLE1<-Weibull_moments(alpha=RMLE1[1],lambda=RMLE1[2])
  
  MoM2<-median_of_means(sortedx,korder=2)
  MoM3<-median_of_means(sortedx,korder=3)
  MoM4<-median_of_means(sortedx,korder=4)
  MoM5<-median_of_means(sortedx,korder=5)
  MoM6<-median_of_means(sortedx,korder=6)
  
  MoRMall<-mHLM_all(x=sortedx,max_dim=6,orderlists=orderlist1_hllarge_asymptotic_rand,boot=TRUE,quasi=FALSE,largesize=largesize)
  
  MoRMall<-unlist(MoRMall)
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[3899],imoments1[3900],imoments1[3901])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)-imoments1[3699])/momentssd[1],
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-imoments1[3773])/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-imoments1[3827])/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-imoments1[3869])/momentssd[4],
    skewbias=abs(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])-imoments1[3903]),
    kurtbias=abs(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902])-imoments1[3902])
  )
  
  allrawmo1<-c(first=(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[(1779):(2498)],imoments1[3699:3736],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoRMall=MoRMall)),
               second=(c(imoments1[2],momentsx[2],imoments1[2499:3018],imoments1[3773:3800],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])),
               third=(c(imoments1[3],momentsx[3],imoments1[3019:3418],imoments1[3827:3848],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])),
               fourth=(c(imoments1[4],momentsx[4],imoments1[3419:3698],imoments1[3869:3884],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])),
               skew=(c(imoments1[5],momentsx[7],imoments1[737:1778],imoments1[3903])),
               kurt=(c(imoments1[6],momentsx[8],imoments1[7:736],imoments1[3902]))
  )
  
  medianmoments<-c(imoments1[3706],imoments1[3780],imoments1[3834],imoments1[3875])
  standardizedm<-c(imoments1[3706]/momentssd[1],imoments1[3780]/momentssd[2],imoments1[3834]/momentssd[3],imoments1[3875]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}
write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_gnorm_raw_Process_max",largesize,".csv", sep = ","), row.names = FALSE)


stopCluster(cl)
registerDoSEQ()

