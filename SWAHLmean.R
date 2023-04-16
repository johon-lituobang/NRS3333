#Copyright 2023 Tuban Lee
#This is a test package that is currently under review in PNAS, please do not share it.


SWHLM<-function(x,max_dim=6,orderlists=NULL,percentage=1/24,batch = "auto", boot = TRUE){
  sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(sortedx)
  
  result_list<-list()
  
  bootstrappedsample1<-bootbatch1(orderlist11=orderlists,sortedx=sortedx,dimension=max_dim)
  if(1 - (1 - percentage)^max_dim>0.5){
    max_dim<-floor(log(1-1/2)/log(1-percentage))
  }
  if(lengthx>40 || boot){
    for(dim in 2:max_dim){
      
      if(dim==2){
        dp_hl<-hlapply1(bootstrappedsample1)
      }else{
        function_name<-paste0("hl", dim, "apply1")
        function_to_call<-get(function_name)
        dp_hl<-function_to_call(bootstrappedsample1)
      }
      
      dp_hl_sorted<-Rfast::Sort(x=dp_hl,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
      
      SWHLM_dim<-SWA(x=dp_hl_sorted,percentage=1-(1-percentage)^dim,sorted=TRUE)
      
      result_list[[paste0("SWHLM",dim)]]<-SWHLM_dim
    }
  }
  
  return(result_list)
}

SWAHLmean<-function (x,orderlist1_sorted2=NULL,orderlist1_sorted3=NULL,orderlist1_sorted4=NULL,orderlist1_sorted5=NULL,orderlist1_sorted6=NULL,batch="auto",boot=TRUE){
  sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

  lengthx<-length(sortedx)

  if(lengthx>40 || boot){
    
    bootstrappedsample2<-bootbatch1(orderlist11 = orderlist1_sorted2,sortedx=sortedx,dimension=2)

    dp2hl<-hlapply1(bootstrappedsample2)
    dp2hlsorted<-Rfast::Sort(x=dp2hl,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

    bootstrappedsample2<-c()

    bootstrappedsample3<-bootbatch1(orderlist11=orderlist1_sorted3,sortedx=sortedx,dimension=3)
    dp3hl3<-hl3apply1(bootstrappedsample3)
    dp3hl3sorted<-Rfast::Sort(x=dp3hl3,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

    bootstrappedsample3<-c()
    
    bootstrappedsample4<-bootbatch1(orderlist11 = orderlist1_sorted4,sortedx=sortedx,dimension=4)
    dp4hl4<-hl4apply1(bootstrappedsample4)
    dp4hl4sorted<-Rfast::Sort(x=dp4hl4,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    bootstrappedsample4<-c()
    
    length1<-min(c(length(dp2hl),length(dp3hl3)))
    dp2hl215<-dp2hl[1:round(length1*(1-((2*log(2)-log(3))/(3*log(2)-log(7)))+floor((2*log(2)-log(3))/(3*log(2)-log(7)))))]
    
    dp2hl315<-dp3hl3[(round(length1*(1-((2*log(2)-log(3))/(3*log(2)-log(7)))+floor((2*log(2)-log(3))/(3*log(2)-log(7)))))+1):length1]
    
    dp1215<-c(dp2hl215,dp2hl315)
    dp1215<-Rfast::Sort(x=dp1215,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    bootstrappedsample5<-bootbatch1(orderlist11 = orderlist1_sorted5,sortedx=sortedx,dimension=5)
    dp4hl5<-hl5apply1(bootstrappedsample5)
    dp4hl5sorted<-Rfast::Sort(x=dp4hl5,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)

    bootstrappedsample5<-c()
    
    bootstrappedsample6<-bootbatch1(orderlist11 = orderlist1_sorted6,sortedx=sortedx,dimension=6)
    dp4hl6<-hl6apply1(bootstrappedsample6)
    dp4hl6sorted<-Rfast::Sort(x=dp4hl6,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    
    length2<-min(c(length(dp4hl5),length(dp4hl6)))
    dp2hl519<-dp4hl5[1:round(length2*(1-((log(2))/(3*log(2)-log(7)))+floor((log(2))/(3*log(2)-log(7)))))]
    
    dp2hl619<-dp4hl6[(round(length2*(1-((log(2))/(3*log(2)-log(7)))+floor((log(2))/(3*log(2)-log(7)))))+1):length2]
    
    dp1519<-c(dp2hl519,dp2hl619)
    dp1519<-Rfast::Sort(x=dp1519,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    bootstrappedsample6<-c()
    
  }

  sq4hl<-c(quantilefunction(x=dp1215,quatiletarget=1/4,sorted=TRUE),quantilefunction(x=dp1215,quatiletarget=3/4,sorted=TRUE))
  mHLM1<-mediansorted(sortedx = dp1519,lengthx = length(dp1519))
  
  SWA_hl6_8<-SWA(x=dp4hl6sorted,percentage=1-(1-1/8)^6,batch=batch,sorted=TRUE,rand=TRUE)
  
  SWA_hl5_8<-SWA(x=dp4hl5sorted,percentage=1-(1-1/8)^5,batch=batch,sorted=TRUE,rand=TRUE)
  
  SWA_hl4_8<-SWA(x=dp4hl4sorted,percentage=1-(1-1/8)^4,batch=batch,sorted=TRUE,rand=TRUE)
  
  SWA_hl3_8<-SWA(x=dp3hl3sorted,percentage=1-(1-1/8)^3,batch=batch,sorted=TRUE,rand=TRUE)
  
  SWA_hl2_8<-SWA(x=dp2hlsorted,percentage=1-(1-1/8)^2,batch=batch,sorted=TRUE,rand=TRUE)
  
  sq4hl2<-c(quantilefunction(x=dp2hlsorted,quatiletarget=1/4,sorted=TRUE),quantilefunction(x=dp2hlsorted,quatiletarget=3/4,sorted=TRUE))
  
  SWA_hl_16<-SWA(x=dp2hlsorted,percentage=1/8,batch=batch,sorted=TRUE,rand=TRUE)

  finallall<-c(MHHLM8=sum(sq4hl)/2,
               mHLM8=mHLM1,MHHLM2=sum(sq4hl2)/2,SWA_hl6_8=SWA_hl6_8,SWA_hl5_8=SWA_hl5_8,
               SWA_hl4_8=SWA_hl4_8,
               SWA_hl3_8=SWA_hl3_8,
               SWA_hl2_8=SWA_hl2_8,SWA_hl_16=SWA_hl_16)
  return(c(finallall))
}


trimmed_mean <- function(x, percentage,sorted=FALSE) {
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  n <- length(x_ordered)
  k <- round(n*percentage)
  c(TM=matrixStats::sum2(x_ordered[(k+1):(n-k)])/(n-2*k))
}
Winsorized_mean <- function(x, percentage,sorted=FALSE) {
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  n <- length(x_ordered)
  k <- round(n*percentage)
  (WM=(matrixStats::sum2(x_ordered[(k+1):(n-k)])+x_ordered[(k+1)]*k+x_ordered[(n-k)]*k)/(n))
}
Block_Winsorized_mean <- function(x, percentage,sorted=FALSE) {
  if(sorted){
    x_ordered<-x
  }else{
    x_ordered<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  n <- length(x_ordered)
  k <- round(n*percentage)
  (BWM=(matrixStats::sum2(x_ordered[(k+1):(n-k)])+matrixStats::sum2(x_ordered[(k+1):(2*k)])+matrixStats::sum2(x_ordered[(n-2*k+1):(n-k)]))/(n))
}
median_of_means <- function(x,korder=2) {
  n_blocks=ceiling(length(x)/korder)
  if (korder <= 1) {
    return(median(x))
  }
  data_mixed <- sample(x)
  index_vec <- rep(1:n_blocks, each = ceiling(length(x) / n_blocks))[1:length(x)]
  data_groups <- split(data_mixed, index_vec)
  block_means <- sapply(data_groups, mean)
  return(c(MoM=median(block_means,na.rm = TRUE)))
}

mHLM<-function(x,orderlist1=NULL,dimension=4,boot=TRUE,quasi=TRUE,largesize=1.8*10^4,interval=16){
  sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(sortedx)
  if(lengthx>40 ||boot){
    if(!quasi&&is.null(orderlist1)){
      unibatchran<-matrix(randtoolbox::SFMT(largesize*3*dimension),ncol=dimension)
      
      orderlist1<-createorderlist(quni1=unibatchran,size=lengthx,interval=interval,dimension=dimension)
      orderlist1<-orderlist1[1:largesize,]
    }else if(quasi&&is.null(orderlist1)){
      quasiunisobol_1<-randtoolbox::sobol(n=largesize*3, dim = dimension, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                      mixed = FALSE, method = "C", start = 1)
      orderlist1<-createorderlist(quni1=quasiunisobol_1,size=lengthx,interval=interval,dimension=dimension)
      orderlist1<-orderlist1[1:largesize,]
    }
    bootstrappedsample1<-bootbatch1(orderlist11=orderlist1,sortedx=sortedx,dimension=dimension)
    
    dphlall<-apply_hl(bootstrappedsample1)
    bootstrappedsample1<-c()
    dphlallsorted<-Rfast::Sort(x=dphlall,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    dphlall<-c()
    finallall<-c(mHLM_boot=mediansorted(sortedx=dphlallsorted,lengthx =largesize ))
  }else{
    combinationsample1<-t(as.data.frame(Rfast::comb_n(sortedx,dimension)))
    combinationsample1<-as.matrix(combinationsample1)
    dphlall<-apply_hl(combinationsample1)
    combinationsample1<-c()
    dphlallsorted<-Rfast::Sort(x=dphlall,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    dphlall<-c()
    combinationsample1<-c()
    finallall<-c(mHLM=mediansorted(sortedx=dphlallsorted,lengthx =length(dphlallsorted)))
  }
  return(c(finallall))
}

mHLM_all<-function(x,max_dim=6,orderlists=NULL,boot=TRUE,quasi=TRUE,largesize=1.8*10^4){
  sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(sortedx)
  
  result_list<-list()
  if(!quasi&&is.null(orderlists)){
    unibatchran<-matrix(randtoolbox::SFMT(largesize*3*dimension),ncol=dimension)
    
    orderlists<-createorderlist(quni1=unibatchran,size=lengthx,interval=interval,dimension=dimension)
    orderlists<-orderlists[1:largesize,]
  }else if(quasi&&is.null(orderlists)){
    quasiunisobol_1<-randtoolbox::sobol(n=largesize*3, dim = dimension, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                        mixed = FALSE, method = "C", start = 1)
    orderlists<-createorderlist(quni1=quasiunisobol_1,size=lengthx,interval=interval,dimension=dimension)
    orderlists<-orderlists[1:largesize,]
  }
  bootstrappedsample1<-bootbatch1(orderlist11=orderlists,sortedx=sortedx,dimension=max_dim)
  if(1 - (1 - percentage)^max_dim>0.5){
    max_dim<-floor(log(1-1/2)/log(1-percentage))
  }
  if(lengthx>40 || boot){
    for(dim in 2:max_dim){
      
      if(dim==2){
        dp_hl<-hlapply1(bootstrappedsample1)
      }else{
        function_name<-paste0("hl", dim, "apply1")
        function_to_call<-get(function_name)
        dp_hl<-function_to_call(bootstrappedsample1)
      }
      
      dp_hl_sorted<-Rfast::Sort(x=dp_hl,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
      
      mHLM_dim<-mediansorted(x=dp_hl_sorted,lengthx=largesize)
      
      result_list[[paste0("mHLM",dim)]]<-mHLM_dim
    }
  }
  
  return(result_list)
}
