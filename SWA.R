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
    if(block_size[i]>0){
      block <- x_ordered[idx:(idx + block_size[i] - 1)] 
    }else{
      block <- 0
    }
    block_sums[i] <- matrixStats::sum2(block)
    idx <- idx + block_size[i]
  }
  return(block_sums)
}

SWA<-function (x,percentage,blocknumber=8,batch="auto",sorted=FALSE,rand=TRUE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(20000/(lengthx))+9
  }
  remove_random_values <- function(row,lengthx,reminder) {
    if(reminder==0){
      return(row)
    }else{
      indices_to_remove <- sample(lengthx,reminder)
      row <- row[-indices_to_remove]
      return(row)
    }
  }
  Ksamples<-lengthx*percentage
  IntKsamples<-floor(Ksamples)
  if(percentage==1/2){
    median1=mediansorted(sortedx=x,lengthx)
    return(c(mean=mean(x),median=median1,median=median1))
  }
  if(percentage>1/6&&percentage<=1/4){
    if (rand&&IntKsamples*(1/percentage)!=lengthx){
      target1<-floor(IntKsamples*(1/percentage))
      batcha<-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batcha!=0){
        batchresults1<-NULL
        for (i in 1:batcha) {
          batchresults1<-rbind(batchresults1,SWAall_mh(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=5,sorted=sorted))
        }
      }else{
        batchresults1<-NULL
      }
      
      target1<-ceiling(IntKsamples*(1/percentage))
      batchb<-batch-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batchb!=0){
        batchresults2<-NULL
        for (i in 1:batchb) {
          batchresults2<-rbind(batchresults2,SWAall_mh(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=5,sorted=sorted))
        }
      }else{
        batchresults2<-NULL
      }
      
      return(matrixStats::colMeans2(rbind(batchresults1,batchresults2),useNames=TRUE))
    }else{
      target1<-lengthx
      Groupmean<-SWAall_mh(x,IntKsamples,target1,blocknumber=5,sorted=sorted)
      return(Groupmean)}
  }else if(percentage>1/4){
    if (rand&&IntKsamples*(1/percentage)!=lengthx){
      target1<-floor(IntKsamples*(1/percentage))
      batcha<-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batcha!=0){
        batchresults1<-NULL
        for (i in 1:batcha) {
          batchresults1<-rbind(batchresults1,SWAall_wt(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=3,sorted=sorted))
        }
      }else{
        batchresults1<-NULL
      }
      
      target1<-ceiling(IntKsamples*(1/percentage))
      batchb<-batch-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batchb!=0){
        batchresults2<-NULL
        for (i in 1:batchb) {
          batchresults2<-rbind(batchresults2,SWAall_wt(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=3,sorted=sorted))
        }
      }else{
        batchresults2<-NULL
      }
      
      return(matrixStats::colMeans2(rbind(batchresults1,batchresults2),useNames=TRUE))
      
    }else{
      target1<-lengthx
      Groupmean<-SWAall_wt(x,IntKsamples,target1,blocknumber=3,sorted=sorted)}
    return(Groupmean)
  }else if(percentage>1/8&&percentage<=1/6){
    if (rand&&IntKsamples*(1/percentage)!=lengthx){
      target1<-floor(IntKsamples*(1/percentage))
      batcha<-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batcha!=0){
        batchresults1<-NULL
        for (i in 1:batcha) {
          batchresults1<-rbind(batchresults1,SWAallbm(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=7,sorted=sorted))
        }
      }else{
        batchresults1<-NULL
      }
      
      target1<-ceiling(IntKsamples*(1/percentage))
      batchb<-batch-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batchb!=0){
        batchresults2<-NULL
        for (i in 1:batchb) {
          batchresults2<-rbind(batchresults2,SWAallbm(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=7,sorted=sorted))
        }
      }else{
        batchresults2<-NULL
      }
      return(matrixStats::colMeans2(rbind(batchresults1,batchresults2),useNames=TRUE))
    }else{
      target1<-lengthx
      Groupmean<-SWAallbm(x,IntKsamples,target1,blocknumber=7,sorted=sorted)
      return(Groupmean)}
  }else if(percentage<=1/8){
    if (rand&&IntKsamples*(1/percentage)!=lengthx){
      target1<-floor(IntKsamples*(1/percentage))
      batcha<-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batcha!=0){
        batchresults1<-NULL
        for (i in 1:batcha) {
          batchresults1<-rbind(batchresults1,SWAall(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=ceiling((target1/IntKsamples)),sorted=sorted))
        }
      }else{
        batchresults1<-NULL
      }
      
      target1<-ceiling(IntKsamples*(1/percentage))
      batchb<-batch-round(((IntKsamples*(1/percentage))%%1)*batch)
      if(batchb!=0){
        batchresults2<-NULL
        for (i in 1:batchb) {
          batchresults2<-rbind(batchresults2,SWAall(remove_random_values(as.numeric(x),lengthx,lengthx-target1),IntKsamples,target1,blocknumber=ceiling((target1/IntKsamples)),sorted=sorted))
        }
      }else{
        batchresults2<-NULL
      }
      return(matrixStats::colMeans2(rbind(batchresults1,batchresults2),useNames=TRUE))
    }else{
      target1<-lengthx
      Groupmean<-SWAall(x,IntKsamples,target1,blocknumber=ceiling((target1/IntKsamples)),sorted=sorted)
      return(Groupmean)}
  }
}

rqm<-function(x,sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  mmm1raw<-SWA8(x=x,interval=8,batch="auto",sorted = TRUE)
  
  mmm1_BM_rm_exp1<-mmmprocessrm(SWA=mmm1raw[2],median=mmm1raw[8],drm=0.37523)
  
  mmm1_BM_qm_exp1<-mmmprocessqm(x=x,percentage=1/8,pc1=CDF(x=x, xevaluated =mmm1raw[2], sorted = TRUE),dqm=0.32128)
  
  c(mmm1_BM_rm_exp1,mmm1_BM_qm_exp1)
}

midhinge<-function(x,sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  sq4<-c(quantilefunction(x=x,quatiletarget=1/4,sorted=TRUE),quantilefunction(x=x,quatiletarget=3/4,sorted=TRUE))
  
  c(midhinge=sum(sq4)/2)
}
stratified_quantile_mean <- function(x,epsilon,gamma,sorted=TRUE) {
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
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
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
  intervalsum1<-intervalsum(x_ordered=x,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
  
  names(intervalsum1)<-NULL
  
  median1<-mediansorted(sortedx=x,lengthx=target1)
  names(median1)<-NULL
  
  TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
  
  winsor<-(c(x[(IntKsamples+1)],x[(((target1-IntKsamples)))]))
  
  winsor1<-sum(winsor)*IntKsamples
  
  WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
  
  mean1<-(sum(intervalsum1))/(target1)
  
  Groupmean<-c(mean=mean1,WM=WM1,TM=TM1,median=median1)
  
  return(Groupmean)
}

SWAall_mh<-function(x,IntKsamples,target1,blocknumber,sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
  intervalsum1<-intervalsum(x_ordered=x,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
  
  names(intervalsum1)<-NULL
  
  median1<-mediansorted(sortedx=x,lengthx=target1)
  names(median1)<-NULL
  
  TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
  
  winsor<-(c(x[(IntKsamples+1)],x[(((target1-IntKsamples)))]))
  
  winsor1<-sum(winsor)*IntKsamples
  
  WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
  
  BWM1<-(TM1*(target1-IntKsamples*2)+matrixStats::sum2(x[(IntKsamples+1):(2*IntKsamples)])+matrixStats::sum2(x[(target1-2*IntKsamples+1):(target1-IntKsamples)]))/(target1)
  
  MH1<-((x[(IntKsamples)+1]+x[target1-(IntKsamples)]))
  if(((target1/IntKsamples)/4)%%1==0){
    SQM1<-MH1/2
  }else{
    SQM1<-(MH1*(IntKsamples*2)+intervalsum1[3])/target1
  }
  
  mean1<-(sum(intervalsum1))/(target1)
  
  Groupmean<-c(mean=mean1,MH=SQM1,WM=WM1,BWM=BWM1,TM=TM1,median=median1)
  
  return(Groupmean)
}

SWAall<-function(x,IntKsamples,target1,blocknumber,sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
  intervalsum1<-intervalsum(x_ordered=x,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
  
  names(intervalsum1)<-NULL
  
  median1<-mediansorted(sortedx=x,lengthx=target1)
  names(median1)<-NULL
  
  TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
  
  winsor<-(c(x[(IntKsamples+1)],x[(((target1-IntKsamples)))]))
  
  winsor1<-sum(winsor)*IntKsamples
  
  WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
  
  BWM1<-(TM1*(target1-IntKsamples*2)+matrixStats::sum2(x[(IntKsamples+1):(2*IntKsamples)])+matrixStats::sum2(x[(target1-2*IntKsamples+1):(target1-IntKsamples)]))/(target1)
  
  middle_blockindex<-ceiling(blocknumber / 2)
  
  m <- floor(1 / (4 * (IntKsamples/target1)))
  result <- 0
  for (i in 1:m) {
    term1 <- x[(IntKsamples)*(2 * i - 1)+1]
    term2 <- x[target1-(IntKsamples)*(2 * i - 1)]
    result <- result + (term1 + term2)
  }
  
  if(((target1/IntKsamples)/4)%%1==0){
    SQM1<-((result)/(m*2))
  }else{
    SQM1<-(((result))*(IntKsamples*2)+matrixStats::sum2(x[((IntKsamples)*(2*(m))+1):(target1-(IntKsamples)*(2*m))]))/target1
  }
  
  BM2result<-0
  
  BM2m<-floor((((target1/IntKsamples)/6)))
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
  BM3m<-floor(target1/IntKsamples/ 8)
  for (i in 1:BM3m) {
    BM3term1 <- intervalsum1[1+(i-1)*4]+intervalsum1[length(intervalsum1)-(1-1+(i-1)*4)]
    BM3term2 <- intervalsum1[1+1+(i-1)*4]+intervalsum1[length(intervalsum1)-(1+(i-1)*4)]
    BM3term3 <- intervalsum1[1+2+(i-1)*4]+intervalsum1[length(intervalsum1)-(1+1+(i-1)*4)]
    BM3term4 <- intervalsum1[1+3+(i-1)*4]+intervalsum1[length(intervalsum1)-(1+2+(i-1)*4)]
    BM3result <- BM3result + (BM3term1*0 + BM3term2*4-BM3term3*2+BM3term4*2)
  }
  if(((target1/IntKsamples)/8)%%1==0){
    BM31<-((BM3result)/(target1))
  }else{
    BM31<-(BM3result+matrixStats::sum2(intervalsum1[(BM3m*4+1):(length(intervalsum1)+1-(BM3m*4+1))]))/target1
  }
  
  mean1<-(sum(intervalsum1))/(target1)
  
  Groupmean<-c(mean=mean1,BM3=BM31,SQM=SQM1,BM2=BM21,WM=WM1,BWM=BWM1,TM=TM1,median=median1)
  
  return(Groupmean)
}

SWAallbm<-function(x,IntKsamples,target1,blocknumber,sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  
  intervalsum1<-intervalsum(x_ordered=x,IntKsamples=IntKsamples,target1=target1,blocknumber=blocknumber)
  
  names(intervalsum1)<-NULL
  
  median1<-mediansorted(sortedx=x,lengthx=target1)
  names(median1)<-NULL
  
  TM1<-matrixStats::sum2(intervalsum1[2:(length(intervalsum1)-1)])/(target1-IntKsamples*2)
  
  winsor<-(c(x[(IntKsamples+1)],x[(((target1-IntKsamples)))]))
  
  winsor1<-sum(winsor)*IntKsamples
  
  WM1<-(TM1*(target1-IntKsamples*2)+winsor1)/(target1)
  
  BWM1<-(TM1*(target1-IntKsamples*2)+matrixStats::sum2(x[(IntKsamples+1):(2*IntKsamples)])+matrixStats::sum2(x[(target1-2*IntKsamples+1):(target1-IntKsamples)]))/(target1)
  
  MH1<-((x[(IntKsamples)+1]+x[target1-(IntKsamples)]))
  
  SQM1<-(MH1*(IntKsamples*2)+matrixStats::sum2(intervalsum1[3:5]))/target1
  
  if (IntKsamples/target1==1/6){
    BM21<-((intervalsum1[2]+intervalsum1[6])/(IntKsamples*2))
  }else{
    BM21<-((intervalsum1[2]+intervalsum1[6])*3+intervalsum1[4])/target1
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
      x<-x
    }else{
      x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    }
    if ((target1/IntKsamples)==9){
      onelist<-(c(x[1:(1*(IntKsamples))],x[(target1-1*IntKsamples+1):(target1)]))
      twolist<-(c(x[(IntKsamples+1):(2*(IntKsamples))],x[(target1-2*IntKsamples+1):(target1-IntKsamples)]))
      threelist<-(c(x[(2*(IntKsamples)+1):(3*(IntKsamples))],x[(6*(IntKsamples)+1):(target1-2*IntKsamples)]))
      fourlist<-(c(x[(3*(IntKsamples)+1):(4*(IntKsamples))],x[(5*(IntKsamples)+1):(target1-3*IntKsamples)]))
      fivelist<-(x[(4*(IntKsamples)+1):(5*(IntKsamples))])
      
      sonelist<-sum(onelist)
      stwolist<-sum(twolist)
      sthreelist<-sum(threelist)
      sfourlist<-sum(fourlist)
      sfivelist<-sum(fivelist)
      
      winsor<-(c(x[(IntKsamples+1)],x[((target1-IntKsamples))]))
      winsor1<-sum(winsor)*IntKsamples
      wm1<-(stwolist+sthreelist+sfourlist+sfivelist+winsor1)/(IntKsamples*(9))
      
      Groupmean<-c(SM9=(sfivelist+stwolist)/(IntKsamples*3),WM9=wm1)
      
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


SWA8<-function (x,interval=8,batch="auto",sorted=FALSE){
  if(sorted){
    x<-x
  }else{
    x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(500000/lengthx)+1
  }
  if(interval!=8){
    return("interval must be 8 ")
  }
  Ksamples<-lengthx/interval
  IntKsamples<-floor(Ksamples)
  target1<-IntKsamples*interval
  smass<-function(x,IntKsamples,target1,sorted=FALSE){
    if(sorted){
      x<-x
    }else{
      x<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    }
    if ((target1/IntKsamples)==8){
      list1<-(c(x[1:(1*(IntKsamples))],x[(7*IntKsamples+1):(target1)]))
      list2<-(c(x[(IntKsamples+1):(2*(IntKsamples))],x[(6*IntKsamples+1):(7*IntKsamples)]))
      list3<-(c(x[(2*(IntKsamples)+1):(3*(IntKsamples))],x[(5*(IntKsamples)+1):(6*IntKsamples)]))
      list4<-c(x[(3*(IntKsamples)+1):(5*(IntKsamples))])
      
      sumlist1<-matrixStats::sum2(list1)
      sumlist2<-matrixStats::sum2(list2)
      sumlist3<-matrixStats::sum2(list3)
      sumlist4<-matrixStats::sum2(list4)
      
      median1<-mediansorted(sortedx=x,lengthx=target1)
      names(median1)<-NULL
      tm8<-(sumlist2+sumlist3+sumlist4)/(IntKsamples*(6))
      
      winsor<-(c(x[(IntKsamples+1)],x[(((target1-IntKsamples)))]))
      
      winsor1<-sum(winsor)*IntKsamples
      
      wm8<-(tm8*IntKsamples*(6)+winsor1)/(IntKsamples*(8))
      
      bwm8<-(tm8*IntKsamples*(6)+sumlist2)/(IntKsamples*(8))
      
      sq8<-c(quantilefunction(x=x,quatiletarget=1/8,sorted=TRUE),quantilefunction(x=x,quatiletarget=3/8,sorted=TRUE),quantilefunction(x=x,quatiletarget=5/8,sorted=TRUE),quantilefunction(x=x,quatiletarget=7/8,sorted=TRUE))
      
      sqm8<-sum(sq8)/4
      
      BM28<-(sumlist2+sumlist3)/(IntKsamples*(4))
      
      BM38<-(sumlist2*4-2*sumlist3+2*sumlist4)/(IntKsamples*(8))
      
      mean1<-(sumlist1+sumlist2+sumlist3+sumlist4)/(IntKsamples*(8))
      
      Groupmean<-c(mean1=mean1,BM38=BM38,sqm8=sqm8,BM28=BM28,wm8=wm8,bwm8=bwm8,tm8=tm8,median1=median1)
      
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


