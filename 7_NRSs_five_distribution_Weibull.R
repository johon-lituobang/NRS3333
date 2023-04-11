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



#bootsize for bootstrap approximation of the distributions of the kernel of U-statistics.
n <- 2048*9*3
(n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (approximating asymptotic)
largesize1<-2048*9

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)

quasiuni<-quasiunisobol

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...

asymptotic_n <- 2048*900*3*2
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (asymptotic bias)
largesize<-2048*900*2

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


orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2_asymptotic,size=largesize,interval=16,dimension=2)
orderlist1_AB2_asymptotic<-orderlist1_AB2[1:largesize,]
orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3_asymptotic,size=largesize,interval=16,dimension=3)
orderlist1_AB3_asymptotic<-orderlist1_AB3[1:largesize,]
orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4_asymptotic,size=largesize,interval=16,dimension=4)
orderlist1_AB4_asymptotic<-orderlist1_AB4[1:largesize,]

quasiuni_sorted2_asymptotic<-c()
quasiuni_sorted3_asymptotic<-c()
quasiuni_sorted4_asymptotic<-c()


asymptotic_n <- 2048*900*3*2
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 10
#large sample size (asymptotic bias)
largesize<-2048*900*2

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_M<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

orderlist1_AB10_quasi_asymptotic<-createorderlist(quni1=quasiuni_M[,1:10],size=largesize,interval=16,dimension=10)
orderlist1_AB10_quasi_asymptotic<-orderlist1_AB10_quasi_asymptotic[1:largesize,]

quasiuni_M<-c()
setSeed(1)

unibatchran_M<-matrix(randtoolbox::SFMT(largesize*3*morder),ncol=morder)

orderlist1_AB10_random_asymptotic<-createorderlist(quni1=unibatchran_M[,1:10],size=largesize,interval=16,dimension=10)
orderlist1_AB10_random_asymptotic<-orderlist1_AB10_random_asymptotic[1:largesize,]


n <- 2048*9*3
(n%%10)==0
# maximum order of moments
morder <- 10
#large sample size (asymptotic bias)
largesize1<-2048*9
samplesize=2048*2
#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_M<-rbind(quasiunisobol)

quasiunisobol<-c()

orderlist1_AB10_quasi<-createorderlist(quni1=quasiuni_M[,1:10],size=samplesize,interval=16,dimension=10)
orderlist1_AB10_quasi<-orderlist1_AB10_quasi[1:largesize1,]

quasiuni_M<-c()
setSeed(1)

orderlist1_AB10_randomall<-c()
for(i in (1:1000)){
  unibatchran_M<-matrix(randtoolbox::SFMT(largesize1*3*morder),ncol=morder)
  
  orderlist1_AB10_random<-createorderlist(quni1=unibatchran_M[,1:10],size=samplesize,interval=16,dimension=10)
  orderlist1_AB10_random<-orderlist1_AB10_random[1:largesize1,]
  
  unibatchran_M<-c()
  
  orderlist1_AB10_randomall<-cbind(orderlist1_AB10_randomall,orderlist1_AB10_random)
}



#load asymptotic d for two parameter distributions
d_values<- read.csv(("d_values.csv"))
I_values<- read.csv(("I_values.csv"))
Imoments_values<- read.csv(("Imoments_values.csv"))


#set the stop criterion
criterionset=1e-10

kurtWeibull<- read.csv(("kurtWeibull_31150.csv"))
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
  SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  imoments1<-unlist(imoments1)
  momentsx<-unbiasedmoments(x=sortedx)
  
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
  MoM7<-median_of_means(sortedx,korder=7)
  MoM8<-median_of_means(sortedx,korder=8)
  MoM9<-median_of_means(sortedx,korder=9)
  MoM10<-median_of_means(sortedx,korder=10)
  MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_random_asymptotic,boot=TRUE,quasi=FALSE,largesize=largesize)
  mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi_asymptotic,boot=TRUE,quasi=TRUE,largesize=largesize)
  
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm)/sqrt(targetvar),
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm)/momentssd[4])
  allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
               second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
               third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
               fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
  )
  
  
  medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
  standardizedm<-c(imoments1[1398]/sqrt(targetvar),imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_Weibull_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtWeibull),c(1:345)]),paste("asymptotic_Weibull",largesize,".csv", sep = ","), row.names = FALSE)

largesize1<-2048*9
samplesize=2048*2
batchsizebase=1000
orderlist1_AB20<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=16,dimension=2)
orderlist1_AB20<-orderlist1_AB20[1:largesize1,]
orderlist1_AB30<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=16,dimension=3)
orderlist1_AB30<-orderlist1_AB30[1:largesize1,]
orderlist1_AB40<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=16,dimension=4)
orderlist1_AB40<-orderlist1_AB40[1:largesize1,]

orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2,size=largesize1,interval=16,dimension=2)
orderlist1_AB2<-orderlist1_AB2[1:largesize1,]
orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3,size=largesize1,interval=16,dimension=3)
orderlist1_AB3<-orderlist1_AB3[1:largesize1,]
orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4,size=largesize1,interval=16,dimension=4)
orderlist1_AB4<-orderlist1_AB4[1:largesize1,]

batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
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
    
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    
    MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_randomall[,(seq(from=1, to=10000,by=10)[batchnumber]):(seq(from=10, to=10000,by=10)[batchnumber])],boot=TRUE,quasi=FALSE,largesize=largesize1)
    mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi,boot=TRUE,quasi=TRUE,largesize=largesize1)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("Weibull_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  SEbatachesmean[]
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,375])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  
  a=allkurtWeibull[batchnumber]
  
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("Weibull_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("Weibull_ABSSE_error.csv", sep = ","), row.names = FALSE)


simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  library(randtoolbox)
  setSeed(1)
  set.seed(1)
  a=allkurtWeibull[60]
  
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,boot=TRUE,quasi=FALSE,largesize=largesize)
    mHLMall<-mHLM10(x=sortedx,boot=TRUE,quasi=TRUE,largesize=largesize)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("finite_Weibull_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[60,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[60,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[60,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[60,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[60,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[60,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[60,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[60,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[60,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[60,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[60,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[60,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[60,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[60,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[60,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[60,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[60,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[60,362]
  
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
  
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_bias_Monte,paste("finite_Weibull_bootstrapsize_raw",samplesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  largesize<-round(batchnumber*2048*9/100)
  a=allkurtWeibull[60]
  
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("finite_Weibull_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[60,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[60,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[60,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[60,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[60,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[60,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[60,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[60,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[60,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[60,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[60,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[60,362]
  
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
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatch_bias_Monte_SE,paste("finite_Weibull_bootstrapsize_raw_error",samplesize,".csv", sep = ","), row.names = FALSE)



kurtgamma<- read.csv(("kurtgamma_31150.csv"))
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
  SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  imoments1<-unlist(imoments1)
  momentsx<-unbiasedmoments(x=sortedx)
  
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
  MoM7<-median_of_means(sortedx,korder=7)
  MoM8<-median_of_means(sortedx,korder=8)
  MoM9<-median_of_means(sortedx,korder=9)
  MoM10<-median_of_means(sortedx,korder=10)
  MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_random_asymptotic,boot=TRUE,quasi=FALSE,largesize=largesize)
  mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi_asymptotic,boot=TRUE,quasi=TRUE,largesize=largesize)
  
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm)/sqrt(targetvar),
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm)/momentssd[4])
  allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
               second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
               third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
               fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
  )
  
  
  medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
  standardizedm<-c(imoments1[1398]/sqrt(targetvar),imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_gamma_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtgamma),c(1:345)]),paste("asymptotic_gamma",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtgamma))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  
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
    
    
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_randomall[,(seq(from=1, to=10000,by=10)[batchnumber]):(seq(from=10, to=10000,by=10)[batchnumber])],boot=TRUE,quasi=FALSE,largesize=largesize1)
    mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi,boot=TRUE,quasi=TRUE,largesize=largesize1)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  
  write.csv(SEbataches,paste("gamma_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  a=allkurtgamma[batchnumber]
  
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("gamma_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("gamma_ABSSE_error.csv", sep = ","), row.names = FALSE)



simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  library(randtoolbox)
  setSeed(1)
  set.seed(1)
  a=allkurtgamma[15]
  
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
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
    
    x<-c(dsgamma(uni=unibatch[,batch1], shape=a/1, rate  = 1))
    
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,boot=TRUE,quasi=FALSE,largesize=largesize)
    mHLMall<-mHLM10(x=sortedx,boot=TRUE,quasi=TRUE,largesize=largesize)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("finite_gamma_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[15,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[15,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[15,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[15,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[15,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[15,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[15,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[15,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[15,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[15,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[15,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[15,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[15,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[15,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[15,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[15,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[15,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[15,362]
  
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
  
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_bias_Monte,paste("finite_gamma_bootstrapsize_raw",samplesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  largesize<-round(batchnumber*2048*9/100)
  a=allkurtgamma[15]
  
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("finite_gamma_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[15,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[15,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[15,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[15,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[15,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[15,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[15,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[15,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[15,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[15,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[15,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[15,362]
  
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
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatch_bias_Monte_SE,paste("finite_gamma_bootstrapsize_raw_error",samplesize,".csv", sep = ","), row.names = FALSE)




kurtPareto<- read.csv(("kurtPareto_91210.csv"))
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
  SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  imoments1<-unlist(imoments1)
  momentsx<-unbiasedmoments(x=sortedx)
  
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
  MoM7<-median_of_means(sortedx,korder=7)
  MoM8<-median_of_means(sortedx,korder=8)
  MoM9<-median_of_means(sortedx,korder=9)
  MoM10<-median_of_means(sortedx,korder=10)
  MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_random_asymptotic,boot=TRUE,quasi=FALSE,largesize=largesize)
  mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi_asymptotic,boot=TRUE,quasi=TRUE,largesize=largesize)
  
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm)/sqrt(targetvar),
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm)/momentssd[4])
  allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
               second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
               third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
               fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
  )
  
  
  medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
  standardizedm<-c(imoments1[1398]/sqrt(targetvar),imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_Pareto_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtPareto),1:345]),paste("asymptotic_Pareto",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtPareto))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  
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
    
    
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_randomall[,(seq(from=1, to=10000,by=10)[batchnumber]):(seq(from=10, to=10000,by=10)[batchnumber])],boot=TRUE,quasi=FALSE,largesize=largesize1)
    mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi,boot=TRUE,quasi=TRUE,largesize=largesize1)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  
  write.csv(SEbataches,paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  
  a=allkurtPareto[batchnumber]
  
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatch_ABSE_SE,paste("Pareto_ABSSE_error.csv", sep = ","), row.names = FALSE)



simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  library(randtoolbox)
  setSeed(1)
  set.seed(1)
  
  a=allkurtPareto[88]
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  
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
    x<-c(dsPareto(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,boot=TRUE,quasi=FALSE,largesize=largesize)
    mHLMall<-mHLM10(x=sortedx,boot=TRUE,quasi=TRUE,largesize=largesize)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("finite_Pareto_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[88,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[88,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[88,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[88,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[88,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[88,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[88,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[88,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[88,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[88,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[88,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[88,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[88,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[88,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[88,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[88,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[88,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[88,362]
  
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
  
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_bias_Monte,paste("finite_Pareto_bootstrapsize_raw",samplesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  largesize<-round(batchnumber*2048*9/100)
  a=allkurtPareto[88]
  
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("finite_Pareto_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[88,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[88,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[88,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[88,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[88,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[88,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[88,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[88,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[88,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[88,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[88,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[88,362]
  
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
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatch_bias_Monte_SE,paste("finite_Pareto_bootstrapsize_raw_error",samplesize,".csv", sep = ","), row.names = FALSE)




kurtlognorm<- read.csv(("kurtlognorm_31150.csv"))
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
  SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  imoments1<-unlist(imoments1)
  momentsx<-unbiasedmoments(x=sortedx)
  
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
  MoM7<-median_of_means(sortedx,korder=7)
  MoM8<-median_of_means(sortedx,korder=8)
  MoM9<-median_of_means(sortedx,korder=9)
  MoM10<-median_of_means(sortedx,korder=10)
  MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_random_asymptotic,boot=TRUE,quasi=FALSE,largesize=largesize)
  mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi_asymptotic,boot=TRUE,quasi=TRUE,largesize=largesize)
  
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm)/sqrt(targetvar),
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm)/momentssd[4])
  allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
               second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
               third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
               fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
  )
  
  
  medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
  standardizedm<-c(imoments1[1398]/sqrt(targetvar),imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}

write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_lognorm_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtlognorm),1:345]),paste("asymptotic_lognorm",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtlognorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
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
    
    
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_randomall[,(seq(from=1, to=10000,by=10)[batchnumber]):(seq(from=10, to=10000,by=10)[batchnumber])],boot=TRUE,quasi=FALSE,largesize=largesize1)
    mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi,boot=TRUE,quasi=TRUE,largesize=largesize1)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  
  write.csv(SEbataches,paste("lognorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  
  a=allkurtlognorm[batchnumber]
  
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("lognorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}



write.csv(simulatedbatch_ABSE_SE,paste("lognorm_ABSSE_error.csv", sep = ","), row.names = FALSE)



simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  library(randtoolbox)
  setSeed(1)
  set.seed(1)
  
  a=allkurtlognorm[59]
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  
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
    x<-c(dslnorm(uni=unibatch[,batch1],  meanlog =0, sdlog = a/1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,boot=TRUE,quasi=FALSE,largesize=largesize)
    mHLMall<-mHLM10(x=sortedx,boot=TRUE,quasi=TRUE,largesize=largesize)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("finite_lognormal_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[59,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[59,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[59,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[59,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[59,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[59,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[59,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[59,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[59,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[59,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[59,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[59,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[59,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[59,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[59,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[59,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[59,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[59,362]
  
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
  
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_bias_Monte,paste("finite_lognormal_bootstrapsize_raw",samplesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  largesize<-round(batchnumber*2048*9/100)
  
  a=allkurtlognorm[59]
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("finite_lognormal_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[59,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[59,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[59,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[59,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[59,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[59,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[59,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[59,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[59,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[59,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[59,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[59,362]
  
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
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatch_bias_Monte_SE,paste("finite_lognormal_bootstrapsize_raw_error",samplesize,".csv", sep = ","), row.names = FALSE)



kurtgnorm<- read.csv(("kurtgnorm_31150.csv"))
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
  SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
  imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB2_asymptotic,orderlist1_sorted30=orderlist1_AB3_asymptotic,orderlist1_sorted40=orderlist1_AB4_asymptotic,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
  imoments1<-unlist(imoments1)
  momentsx<-unbiasedmoments(x=sortedx)
  
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
  MoM7<-median_of_means(sortedx,korder=7)
  MoM8<-median_of_means(sortedx,korder=8)
  MoM9<-median_of_means(sortedx,korder=9)
  MoM10<-median_of_means(sortedx,korder=10)
  MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_random_asymptotic,boot=TRUE,quasi=FALSE,largesize=largesize)
  mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi_asymptotic,boot=TRUE,quasi=TRUE,largesize=largesize)
  
  sortedx<-c()
  momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
  
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm)/sqrt(targetvar),
    secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar)/momentssd[2],
    thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm)/momentssd[3],
    fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm)/momentssd[4])
  allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
               second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
               third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
               fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
  )
  
  
  medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
  standardizedm<-c(imoments1[1398]/sqrt(targetvar),imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
  all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
}
write.csv(simulatedbatch_asymptoticbias,paste("asymptotic_gnorm_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatch_asymptoticbias[1:length(allkurtgnorm),1:345]),paste("asymptotic_gnorm",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtgnorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  
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
    
    
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_randomall[,(seq(from=1, to=10000,by=10)[batchnumber]):(seq(from=10, to=10000,by=10)[batchnumber])],boot=TRUE,quasi=FALSE,largesize=largesize1)
    mHLMall<-mHLM10(x=sortedx,orderlist1=orderlist1_AB10_quasi,boot=TRUE,quasi=TRUE,largesize=largesize1)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  
  write.csv(SEbataches,paste("gnorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[batchnumber,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[batchnumber,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[batchnumber,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  library(randtoolbox)
  library(NRSReview)
  setSeed(1)
  set.seed(1)
  
  
  a=allkurtgnorm[batchnumber]
  
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("gnorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[batchnumber,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[batchnumber,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[batchnumber,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[batchnumber,362]
  
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
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatch_ABSE_SE,paste("gnorm_ABSSE_error.csv", sep = ","), row.names = FALSE)



simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  library(randtoolbox)
  setSeed(1)
  set.seed(1)
  
  a=allkurtgnorm[1]
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  
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
    x<-c(dsgnorm(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch=1,sorted=TRUE)
    imoments1<-imoments(x=sortedx,dtype1=1,releaseall=TRUE,standist_d=d_values,standist_I=I_values,standist_Imoments=Imoments_values,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",stepsize=1000,criterion=criterionset,boot=TRUE)
    imoments1<-unlist(imoments1)
    momentsx<-unbiasedmoments(x=sortedx)
    
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
    MoM7<-median_of_means(sortedx,korder=7)
    MoM8<-median_of_means(sortedx,korder=8)
    MoM9<-median_of_means(sortedx,korder=9)
    MoM10<-median_of_means(sortedx,korder=10)
    MoRMall<-mHLM10(x=sortedx,boot=TRUE,quasi=FALSE,largesize=largesize)
    mHLMall<-mHLM10(x=sortedx,boot=TRUE,quasi=TRUE,largesize=largesize)
    
    sortedx<-c()
    momentssd<-c(sd=sqrt(momentsx[2]),imoments1[1522],imoments1[1523],imoments1[1524])
    
    allrawmoBias<-c(
      firstbias=abs(c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall)-targetm),
      secondbias=abs(c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2])-targetvar),
      thirdbias=abs(c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3])-targettm),
      fourbias=abs(c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])-targetfm))
    allrawmo1<-c(first=c(Huberx,SMWM9,imoments1[1],momentsx[1],imoments1[7:52],imoments1[121:166],imoments1[1391:1418],mean_medianMAD1=moments_medianMAD1[1],mean_QE1=moments_QE1[1],mean_RMLE1=moments_RMLE1[1],MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM6=MoM6,MoM7=MoM7,MoM8=MoM8,MoM9=MoM9,MoM10=MoM10,MoRMall=MoRMall,mHLMall=mHLMall),
                 second=c(imoments1[2],momentsx[2],imoments1[53:86],imoments1[167:200],imoments1[1442:1463],var_medianMAD1=moments_medianMAD1[2],var_QE1=moments_QE1[2],var_RMLE1=moments_RMLE1[2]),
                 third=c(imoments1[3],momentsx[3],imoments1[87:106],imoments1[201:220],imoments1[1481:1494],tm_medianMAD1=moments_medianMAD1[3],tm_QE1=moments_QE1[3],tm_RMLE1=moments_RMLE1[3]),
                 fourth=c(imoments1[4],momentsx[4],imoments1[107:120],imoments1[221:234],imoments1[1505:1514],fm_medianMAD1=moments_medianMAD1[4],fm_QE1=moments_QE1[4],fm_RMLE1=moments_RMLE1[4])
    )
    
    
    medianmoments<-c(imoments1[1398],imoments1[1449],imoments1[1486],imoments1[1510])
    standardizedm<-c(imoments1[1398]/momentssd[1],imoments1[1449]/momentssd[2],imoments1[1486]/momentssd[3],imoments1[1510]/momentssd[4])
    all1<-(c(kurtx=kurtx,skewx=skewx,momentsx,allrawmoBias,momentssd,medianmoments,standardizedm=standardizedm,allrawmo1,Huberx,SMWM9,imoments1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("finite_gnorm_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,7:161])^2))/sqrt(targetvar)
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,162:256])^2))/simulatedbatch_asymptoticbias[1,360]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,257:315])^2))/simulatedbatch_asymptoticbias[1,361]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,316:358])^2))/simulatedbatch_asymptoticbias[1,362]
  
  AB1_mean<-abs(colMeans((SEbataches[,7:161])))/sqrt(targetvar)
  
  AB1_var<-abs(colMeans((SEbataches[,162:256])))/simulatedbatch_asymptoticbias[1,360]
  
  AB1_tm<-abs(colMeans((SEbataches[,257:315])))/simulatedbatch_asymptoticbias[1,361]
  
  AB1_fm<-abs(colMeans((SEbataches[,316:358])))/simulatedbatch_asymptoticbias[1,362]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[1,360]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[1,361]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[1,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[1,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[1,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[1,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[1,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[1,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[1,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[1,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[1,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[1,362]
  
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
  
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatch_bias_Monte,paste("finite_gnorm_bootstrapsize_raw",samplesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:100)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)
  largesize<-round(batchnumber*2048*9/100)
  a=allkurtgnorm[1]
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("finite_gnorm_bootstrapsize_raw",batchnumber,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,350])
  
  samplemean_SE1<-samplemeansd_unscaled1/sqrt(targetvar)
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,527])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatch_asymptoticbias[1,360]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,622])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatch_asymptoticbias[1,361]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,681])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatch_asymptoticbias[1,362]
  
  ratiosamplemean1<-c(SEbatachesmean[350])/SEbatachesmean[468]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,350])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/sqrt(targetvar)
  
  ratiosamplevar1<-c(SEbatachesmean[527])/SEbatachesmean[596]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,527])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatch_asymptoticbias[1,360]
  
  ratiosampletm1<-c(SEbatachesmean[622])/SEbatachesmean[663]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,622])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatch_asymptoticbias[1,361]
  
  ratiosamplefm1<-c(SEbatachesmean[681])/SEbatachesmean[710]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,681])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  
  samplefm_SSE1<-samplefmsd1/simulatedbatch_asymptoticbias[1,362]
  
  ratiomean1<-c(SEbatachesmean[371:525])/SEbatachesmean[468]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,371:525]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,371:525])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/sqrt(targetvar)
  
  ratiovar1<-SEbatachesmean[526:620]/SEbatachesmean[596]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,526:620]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatch_asymptoticbias[1,360]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,526:620])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatch_asymptoticbias[1,360]
  ratiotm1<-SEbatachesmean[621:679]/SEbatachesmean[663]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,621:679]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatch_asymptoticbias[1,361]
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,621:679])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatch_asymptoticbias[1,361]
  
  ratiofm1<-SEbatachesmean[680:722]/SEbatachesmean[710]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,680:722]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatch_asymptoticbias[1,362]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,680:722])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatch_asymptoticbias[1,362]
  
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
  allErrors<-c(bootstrapsize=largesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatch_bias_Monte_SE,paste("finite_gnorm_bootstrapsize_raw_error",samplesize,".csv", sep = ","), row.names = FALSE)



stopCluster(cl)
registerDoSEQ()