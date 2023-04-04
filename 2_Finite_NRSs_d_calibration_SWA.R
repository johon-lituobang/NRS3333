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

#bootsize for bootstrap approximation of the distributions of the kernal of U-statistics.
n <- 2048*9*3
(n%%10)==0
# maximum order of moments
morder <- 4
largesize<-2048*9
#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)

quasiuni<-quasiunisobol

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...


#the function, createorderlist, can transform the Sobol sequence into non-repeated integer sequences.

samplesize=2048*2

#the samplesize is the maximum value of this integer sequence
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
quasiuni_sorted2<-c()
quasiuni_sorted3<-c()
quasiuni_sorted4<-c()


kurtWeibull<- read.csv(("kurtWeibull_28260.csv"))
allkurtWeibull<-unlist(kurtWeibull)

batchsizebase=1000

batchsize=(batchsizebase)

n <- samplesize

setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

#Then, start the Monte Simulation
simulatedbatch_bias_Monte<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
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
    x<-c()
    dall<-compalldmoments(x=sortedx,targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm,orderlist1_sorted20=orderlist1_AB20,orderlist1_sorted30=orderlist1_AB30,orderlist1_sorted40=orderlist1_AB40,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,percentage=1/16,batch="auto",boot=TRUE)
    
    sortedx<-c()
    all1<-t(c(samplesize,type=1,kurtx,skewx,dall))
    SEbataches<-rbind(SEbataches,all1)
  }

  write.csv(SEbataches,paste("finite_Weibull_dcalibration_raw_SWA",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)

  SEbataches
}
resultcolnames1<- read.csv(paste("finite_Weibull_dcalibration_raw_SWA",samplesize,9,".csv", sep = ","))

colnames(simulatedbatch_bias_Monte)<-colnames(resultcolnames1)

Monte_Two<-simulatedbatch_bias_Monte[,c(1,2,3,4,347:ncol(simulatedbatch_bias_Monte))]

simulatedbatch_bias_Monte<-c()

write.csv(Monte_Two,paste("finite_Weibull_dcalibration_raw_SWA",samplesize,".csv", sep = ","), row.names = FALSE)

meanall<-foreach(estimators = c(1:ncol(Monte_Two)), .combine = 'cbind') %dopar% {
  Monte_sample_Estimator1 <- data.frame(Size = ( Monte_Two[,1]),
                                        Kurtosis = (Monte_Two[,3]),
                                        Estimator = Monte_Two[,estimators])

  rownames(Monte_sample_Estimator1)<-c()

  Finite_1<-tapply(Monte_sample_Estimator1$Estimator, Monte_sample_Estimator1$Kurtosis, mean)

  Finite_1
}
colnames(meanall)<-colnames(Monte_Two)
write.csv(meanall,paste("finite_Weibull_dcalibration_raw_SWA_mean",samplesize,".csv", sep = ","), row.names = FALSE)

meanall<-as.data.frame(meanall)
musall1<-meanall[,5:ncol(meanall)]
compalldmoments1<-data.frame(matrix(nrow = nrow(meanall), ncol = 0))
for (estimators in (1:57)){
  mmm1<-musall1[,(9*estimators-8):(9*estimators)]
  mmall<-c()
  for (allestimators in (1:nrow(meanall))){
    mmm1one<-as.numeric(mmm1[allestimators,])
    finddmmm1one1<-compalld(mmm1one)
    mmall<-rbind(mmall,finddmmm1one1)
  }
  compalldmoments1<-cbind(compalldmoments1,mmall)
}
Monte_d_Two<-cbind(meanall[,1:4],compalldmoments1[,c(seq(from=1, to=342, by=3))])

Label_Weibull1<- read.csv(("d_label.csv"))

colnames(Monte_d_Two)<-colnames(Label_Weibull1)

write.csv(Monte_d_Two,paste("finite_d_Weibull_SWA",samplesize,".csv", sep = ","), row.names = FALSE)

simulatedbatch_bias_Monte_SE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(NRSReview)

  a=allkurtWeibull[batchnumber]

  n <- samplesize

  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))

  batchsize=batchsizebase

  SEbataches<- read.csv(paste("finite_Weibull_dcalibration_raw_SWA",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  standarderrors1<-apply(SEbataches,2,se_mean)
  standarderrors1
}
Allstandarderror_each<-cbind(meanall[,1:4],simulatedbatch_bias_Monte_SE[,5:ncol(simulatedbatch_bias_Monte_SE)])

All_each1<-meanall[,5:ncol(meanall)]
Allstandarderror_each1<-Allstandarderror_each[,347:ncol(Allstandarderror_each)]

compdall_error<-c()
for (estimators in (1:57)){
  mmm1_error<-Allstandarderror_each1[,(9*estimators-8):(9*estimators)]
  mmm1_11<-All_each1[,(9*estimators-8):(9*estimators)]

  mmall_error<-c()
  for (allestimators in (1:nrow(Allstandarderror_each1))){
    mmm1one_error<-as.numeric(mmm1_error[allestimators,])
    mmm1one<-as.numeric(mmm1_11[allestimators,])

    compdmmm1one1_error<-compderror(expectanalytic_error=mmm1one_error[1],expecttrue_error=mmm1one_error[2],expectboot_error=mmm1one_error[3],SWA1_error=mmm1one_error[4],median1_error=mmm1one_error[5],mx1_error=mmm1one_error[6],quatileexpectanalytic_error=mmm1one_error[7],
                                    quatileexpecttrue_error=mmm1one_error[8],quatileexpectboot_error=mmm1one_error[9],expectanalytic=mmm1one[1],expecttrue=mmm1one[2],expectboot=mmm1one[3],SWA1=mmm1one[4],median1=mmm1one[5],mx1=mmm1one[6],quatileexpectanalytic=mmm1one[7],
                                    quatileexpecttrue=mmm1one[8],quatileexpectboot=mmm1one[9])
    mmall_error<-rbind(mmall_error,compdmmm1one1_error)
  }
  compdall_error<-cbind(compdall_error,mmall_error)
}
Monte_d_Two_error<-cbind(Allstandarderror_each[,1:4],compdall_error[,c(seq(from=1, to=342, by=3))])

colnames(Monte_d_Two_error)<-colnames(Label_Weibull1)

write.csv(Monte_d_Two_error,paste("finite_d_Weibull_SWA_error.csv", sep = ","), row.names = FALSE)


write.csv(Monte_d_Two,paste("finite_d_SWA.csv", sep = ","), row.names = FALSE)

asymptotic_d_Weibull<- read.csv(paste("asymptotic_d_Weibull_SWA.csv", sep = ","))

asymptotic_d_Weibull<- cbind(dtype=rep(1,nrow(asymptotic_d_Weibull)),asymptotic_d_Weibull)
colnames(asymptotic_d_Weibull)<-NULL
asymptotic_d_Weibull<- cbind(Size=rep(2048*900,nrow(asymptotic_d_Weibull)),asymptotic_d_Weibull)

Label_Weibull1<- read.csv(("d_label.csv"))
colnames(asymptotic_d_Weibull)<-colnames(Label_Weibull1)
d_Merged<- rbind(asymptotic_d_Weibull,Monte_d_Two)
Label_Weibull1<- read.csv(("d_label.csv"))
colnames(d_Merged)<-colnames(Label_Weibull1)

write.csv(d_Merged,paste("d_SWA.csv", sep = ","), row.names = FALSE)

registerDoSEQ()

