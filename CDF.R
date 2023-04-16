#Copyright 2023 Tuban Lee
#This is a test package that is currently under review in PNAS, please do not share it.


CDF<-function(x,xevaluated,sorted=FALSE){
  if(!sorted){
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  matrixStats::mean2(x<=xevaluated)
}
