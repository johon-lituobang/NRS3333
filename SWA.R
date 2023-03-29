#Copyright 2023 Tuban Lee
#This is a test package that is currently under review in PNAS, please do not share it.

intervalsum <- function(x_ordered, IntKsamples, target1, blocknumber = 8) {
  total_blocks <- (target1 / IntKsamples)
  block_size <- rep(IntKsamples, blocknumber)
  
  if (total_blocks != blocknumber) {
    middle_block <- ceiling(blocknumber / 2)
    block_size[middle_block] <- target1 - matrixStats::sum2(block_size[-middle_block])
  }
  
  block_sums <- numeric(blocknumber)
  idx <- 1
  for (i in seq_along(block_size)) {
    block <- x_ordered[idx:(idx + block_size[i] - 1)]
    block_sums[i] <- matrixStats::sum2(block)
    idx <- idx + block_size[i]
  }
  
  return(block_sums)
}

SWA<-function (x,percentage,blocknumber=8,batch="auto",sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(500000/lengthx)+1
  }
  
  Ksamples<-lengthx*percentage
  IntKsamples<-floor(Ksamples)
  target1<-round(IntKsamples*(1/percentage))
  if(percentage>1/6&&percentage<=1/4){
    if (Ksamples%%1!=0 ){
      xt<-t(as.data.frame(x))
      xmatrix<-as.data.frame(lapply(xt, rep, batch))
      allmatrix<-t(apply(xmatrix,1,sample,size=target1,replace=FALSE))
      batchresults<-apply(allmatrix,1,SWAall_mh,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber,sorted=FALSE)
      return(matrixStats::rowMeans2(batchresults,useNames=TRUE))
    }
    else{
      Groupmean<-SWAall_mh(x,IntKsamples,target1,blocknumber,sorted=sorted)
      return(Groupmean)}
  }else if(percentage>1/4){
    if (Ksamples%%1!=0 ){
      xt<-t(as.data.frame(x))
      xmatrix<-as.data.frame(lapply(xt, rep, batch))
      allmatrix<-t(apply(xmatrix,1,sample,size=target1,replace=FALSE))
      batchresults<-apply(allmatrix,1,SWAall_wt,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber,sorted=FALSE)
      return(matrixStats::rowMeans2(batchresults,useNames=TRUE))
    }
    else{
      Groupmean<-SWAall_wt(x,IntKsamples,target1,blocknumber,sorted=sorted)
      return(Groupmean)}
  }else if(percentage>1/8&&percentage<=1/6){
    if (Ksamples%%1!=0 ){
      xt<-t(as.data.frame(x))
      xmatrix<-as.data.frame(lapply(xt, rep, batch))
      allmatrix<-t(apply(xmatrix,1,sample,size=target1,replace=FALSE))
      batchresults<-apply(allmatrix,1,SWAallbm,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber,sorted=FALSE)
      return(matrixStats::rowMeans2(batchresults,useNames=TRUE))
    }
    else{
      Groupmean<-SWAallbm(x,IntKsamples,target1,blocknumber,sorted=sorted)
      return(Groupmean)}
  }else if(percentage<=1/8){
    if (Ksamples%%1!=0 ){
      xt<-t(as.data.frame(x))
      xmatrix<-as.data.frame(lapply(xt, rep, batch))
      allmatrix<-t(apply(xmatrix,1,sample,size=target1,replace=FALSE))
      batchresults<-apply(allmatrix,1,SWAall,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber,sorted=FALSE)
      return(matrixStats::rowMeans2(batchresults,useNames=TRUE))
    }
    else{
      Groupmean<-SWAall(x,IntKsamples,target1,blocknumber,sorted=sorted)
      return(Groupmean)}
  }
  
}

stratified_quantile_mean <- function(x,epsilon,gamma,sorted=TRUE) {
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  n <- length(x)
  m <- floor(1 / (4 * epsilon))
  result <- 0
  
  for (i in 1:m) {
    term1 <- quantilefunction(x, quatiletarget = (2 * i - 1) * gamma * epsilon,sorted=TRUE)
    term2 <- quantilefunction(x, quatiletarget = 1 - (2 * i - 1) * epsilon,sorted=TRUE)
    result <- result + 0.5 * (term1 + term2)
  }
  
  return(4 * epsilon * result)
}

SWAall_wt<-function(x,IntKsamples,target1,blocknumber,sorted=FALSE){
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
  intervalsum1<-intervalsum(x_ordered=x_ordered,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
  
  names(intervalsum1)<-NULL
  
  median1<-mediansorted(sortedx=x_ordered,lengthx=target1)
  names(median1)<-NULL
  
  TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
  
  winsor<-(c(x_ordered[(IntKsamples+1)],x_ordered[(((target1-IntKsamples)))]))
  
  winsor1<-sum(winsor)*IntKsamples
  
  WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
  
  mean1<-(sum(intervalsum1))/(target1)
  
  Groupmean<-c(mean=mean1,WM=WM1,TM=TM1,median=median1)
  
  return(Groupmean)
}

SWAall_mh<-function(x,IntKsamples,target1,blocknumber,sorted=FALSE){
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
  intervalsum1<-intervalsum(x_ordered=x_ordered,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
  
  names(intervalsum1)<-NULL
  
  median1<-mediansorted(sortedx=x_ordered,lengthx=target1)
  names(median1)<-NULL
  
  TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
  
  winsor<-(c(x_ordered[(IntKsamples+1)],x_ordered[(((target1-IntKsamples)))]))
  
  winsor1<-sum(winsor)*IntKsamples
  
  WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
  
  BWM1<-(TM1*(target1-IntKsamples*2)+matrixStats::sum2(x_ordered[(IntKsamples+1):(2*IntKsamples)])+matrixStats::sum2(x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))/(target1)
  
  middle_blockindex<-ceiling(blocknumber / 2)
  
  m <- floor(1 / (4 * (IntKsamples/target1)))
  result <- 0
  for (i in 1:m) {
    term1 <- x_ordered[(IntKsamples)*(2 * i - 1)+1]
    term2 <- x_ordered[target1-(IntKsamples)*(2 * i - 1)]
    result <- result + (term1 + term2)
  }
  
  if(((target1/IntKsamples)/4)%%1==0){
    SQM1<-((result)/(m*2))
  }else{
    SQM1<-(((result)/(m*2))*(IntKsamples*(blocknumber-1))+intervalsum1[middle_blockindex])/target1
  }
  
  mean1<-(sum(intervalsum1))/(target1)
  
  Groupmean<-c(mean=mean1,MH=SQM1,WM=WM1,BWM=BWM1,TM=TM1,median=median1)
  
  return(Groupmean)
}
SWAall<-function(x,IntKsamples,target1,blocknumber,sorted=FALSE){
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
  intervalsum1<-intervalsum(x_ordered=x_ordered,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
  
  names(intervalsum1)<-NULL
  
  median1<-mediansorted(sortedx=x_ordered,lengthx=target1)
  names(median1)<-NULL
  
  TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
  
  winsor<-(c(x_ordered[(IntKsamples+1)],x_ordered[(((target1-IntKsamples)))]))
  
  winsor1<-sum(winsor)*IntKsamples
  
  WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
  
  BWM1<-(TM1*(target1-IntKsamples*2)+matrixStats::sum2(x_ordered[(IntKsamples+1):(2*IntKsamples)])+matrixStats::sum2(x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))/(target1)
  
  middle_blockindex<-ceiling(blocknumber / 2)
  
  m <- floor(1 / (4 * (IntKsamples/target1)))
  result <- 0
  for (i in 1:m) {
    term1 <- x_ordered[(IntKsamples)*(2 * i - 1)+1]
    term2 <- x_ordered[target1-(IntKsamples)*(2 * i - 1)]
    result <- result + (term1 + term2)
  }
  
  if(((target1/IntKsamples)/8)%%1==0){
    SQM1<-((result)/(m*2))
  }else{
    SQM1<-(((result)/(m*2))*(IntKsamples*(blocknumber-1))+intervalsum1[middle_blockindex])/target1
  }
  
  BM2result<-0
  
  BM2m<-floor(((floor(target1/IntKsamples)/2))/3)
  for (i in 1:BM2m) {
    term1 <- intervalsum1[(3 * i - 1)]
    term2 <- intervalsum1[length(intervalsum1)+1-(3 * i - 1)]
    BM2result <- BM2result + (term1 + term2)
  }
  
  if ((floor(target1/IntKsamples)/2)%%3==0){
    BM21<-((BM2result)/(IntKsamples*BM2m*2))
  }else{
    BM21<-((BM2result)*3+matrixStats::sum2(intervalsum1[(BM2m*3+1):(length(intervalsum1)+1-(BM2m*3+1))]))/target1
  }
  
  BM3result<-0
  
  for (i in 1:floor(blocknumber / 8)) {
    BM3term1 <- intervalsum1[1+(i-1)*4]+intervalsum1[length(intervalsum1)-(1-1+(i-1)*4)]
    BM3term2 <- intervalsum1[1+1+(i-1)*4]+intervalsum1[length(intervalsum1)-(1+(i-1)*4)]
    BM3term3 <- intervalsum1[1+2+(i-1)*4]+intervalsum1[length(intervalsum1)-(1+1+(i-1)*4)]
    BM3term4 <- intervalsum1[1+3+(i-1)*4]+intervalsum1[length(intervalsum1)-(1+2+(i-1)*4)]
    BM3result <- BM3result + (BM3term1*0 + BM3term2*4-BM3term3*2+BM3term4*2)
  }
  if(((target1/IntKsamples)/8)%%1==0){
    BM31<-((BM3result)/(target1))
  }else{
    BM31<-(BM3result+intervalsum1[middle_blockindex])/target1
  }
  
  mean1<-(sum(intervalsum1))/(target1)
  
  Groupmean<-c(mean=mean1,BM3=BM31,SQM=SQM1,BM2=BM21,WM=WM1,BWM=BWM1,TM=TM1,median=median1)
  
  return(Groupmean)
}

SWAallbm<-function(x,IntKsamples,target1,blocknumber,sorted=FALSE){
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
    intervalsum1<-intervalsum(x_ordered=x_ordered,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
    
    names(intervalsum1)<-NULL
    
    median1<-mediansorted(sortedx=x_ordered,lengthx=target1)
    names(median1)<-NULL
    
    TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
    
    winsor<-(c(x_ordered[(IntKsamples+1)],x_ordered[(((target1-IntKsamples)))]))
    
    winsor1<-sum(winsor)*IntKsamples
    
    WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
    
    BWM1<-(TM1*(target1-IntKsamples*2)+matrixStats::sum2(x_ordered[(IntKsamples+1):(2*IntKsamples)])+matrixStats::sum2(x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))/(target1)
    
    middle_blockindex<-ceiling(blocknumber / 2)
    
    m <- floor(1 / (4 * (IntKsamples/target1)))
    result <- 0
    for (i in 1:m) {
      term1 <- x_ordered[(IntKsamples)*(2 * i - 1)+1]
      term2 <- x_ordered[target1-(IntKsamples)*(2 * i - 1)]
      result <- result + (term1 + term2)
    }
    
    if(((target1/IntKsamples)/8)%%1==0){
      SQM1<-((result)/(m*2))
    }else{
      SQM1<-(((result)/(m*2))*(IntKsamples*(blocknumber-1))+intervalsum1[middle_blockindex])/target1
    }
    
    BM2result<-0
    
    BM2m<-floor(((floor(target1/IntKsamples)/2))/3)
    for (i in 1:BM2m) {
      term1 <- intervalsum1[(3 * i - 1)]
      term2 <- intervalsum1[length(intervalsum1)+1-(3 * i - 1)]
      BM2result <- BM2result + (term1 + term2)
    }
    
    if ((floor(target1/IntKsamples)/2)%%3==0){
      BM21<-((BM2result)/(IntKsamples*BM2m*2))
    }else{
      BM21<-((BM2result)*3+matrixStats::sum2(intervalsum1[(BM2m*3+1):(length(intervalsum1)+1-(BM2m*3+1))]))/target1
    }
    
    
    mean1<-(sum(intervalsum1))/(target1)
    
    Groupmean<-c(mean=mean1,SQM=SQM1,BM2=BM21,WM=WM1,BWM=BWM1,TM=TM1,median=median1)
    
    return(Groupmean)
}

SWA9<-function (x,interval=9,batch="auto",sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(500000/lengthx)+1
  }
  if(interval!=9){
    return("interval must be 9 ")
  }
  Ksamples<-lengthx/interval
  IntKsamples<-floor(Ksamples)
  target1<-IntKsamples*interval
  smass<-function(x,IntKsamples,target1,sorted=FALSE){
    if(sorted){
      x_ordered<-x
    }else{
      x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    }
    if ((target1/IntKsamples)==9){
      onelist<-(c(x_ordered[1:(1*(IntKsamples))],x_ordered[(target1-1*IntKsamples+1):(target1)]))
      twolist<-(c(x_ordered[(IntKsamples+1):(2*(IntKsamples))],x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))
      threelist<-(c(x_ordered[(2*(IntKsamples)+1):(3*(IntKsamples))],x_ordered[(6*(IntKsamples)+1):(target1-2*IntKsamples)]))
      fourlist<-(c(x_ordered[(3*(IntKsamples)+1):(4*(IntKsamples))],x_ordered[(5*(IntKsamples)+1):(target1-3*IntKsamples)]))
      fivelist<-(x_ordered[(4*(IntKsamples)+1):(5*(IntKsamples))])
      
      sonelist<-sum(onelist)
      stwolist<-sum(twolist)
      sthreelist<-sum(threelist)
      sfourlist<-sum(fourlist)
      sfivelist<-sum(fivelist)
      
      winsor<-(c(x_ordered[(IntKsamples+1)],x_ordered[((target1-IntKsamples))]))
      winsor1<-sum(winsor)*IntKsamples
      wm1<-(stwolist+sthreelist+sfourlist+sfivelist+winsor1)/(IntKsamples*(9))
      
      Groupmean<-c(sm=(sfivelist+stwolist)/(IntKsamples*3),wm=wm1)
      
      return(Groupmean)
    }
    else{
      return("Not supported yet.")
    }
  }
  if (Ksamples%%1!=0 ){
    xt<-t(as.data.frame(x))
    xmatrix<-as.data.frame(lapply(xt, rep, batch))
    allmatrix<-t(apply(xmatrix,1,sample,size=target1,replace=FALSE))
    batchresults<-apply(allmatrix,1,smass,IntKsamples=IntKsamples,target1=target1,sorted=FALSE)
    return(matrixStats::rowMeans2(batchresults,useNames=TRUE))
  }
  else{
    Groupmean<-smass(x,IntKsamples,target1,sorted=sorted)
    return(Groupmean)}
}




x<-rexp(20480)
SWA (x,percentage=1/16,blocknumber=16,batch="auto",sorted=FALSE)
x<-rexp(20480)
SWA (x,percentage=1/8,blocknumber=8,batch="auto",sorted=FALSE)
31/256

x<-rexp(20480)
SWA (x,percentage=31/256,blocknumber=9,batch="auto",sorted=FALSE)

x<-rexp(20480)
SWA (x,percentage=721/4096,blocknumber=5,batch="auto",sorted=FALSE)

x<-rexp(20480)
SWA (x,percentage=14911/65536,blocknumber=5,batch=2,sorted=FALSE)


x<-rexp(20480)
SWA (x,percentage=5386591/16777216,blocknumber=3,batch=2,sorted=FALSE)


x<-rexp(20480)
SWA (x,percentage=1732076671/4294967296,blocknumber=3,batch=2,sorted=FALSE)

x<-rexp(20480)
SWA (x,percentage=30276117361/68719476736,blocknumber=3,batch=2,sorted=FALSE)

