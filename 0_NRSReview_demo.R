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

#this is a demo file explaining this package, NRS.

#This file is totally not relevant to the reviewing of NRS. I will gradually update one I have time..

#you can simulate a distribution with any sample size (better larger than 1000, since the finite sample bias hasn't not been integrated into this package yet.
x<-rexp(4307)

#all robust location estimators are highly bias, e.g.

#SWA function will return all six SWAs with breakdown point 1/8 plus mean and median.

SWA(x)

#compared to recombined mean and quantile mean
rqm(x)

#the principle can be extended to central moments
SWAmoments(x)

#invariant moments
imoments(x)

#invariant standardized moments

istandardizedmoments(x)

#compared to standardized moments
standardizedmoments(x)

#you should ensure that the sample you provided have kurtosis smaller than 26, this is the current limit, but will improved soon.

#finally, median standardized moments are also provided, that is very helpful for studying distributions with infinite moments.
medianstandardizedmoments(x)

#for example
x<-rPareto(21307,shape=1, scale=1)

medianstandardizedmoments(x)
#sample moments are invalid.
standardizedmoments(x)

#median standardized moments are very good at estimating the parameter changes.
x<-rPareto(21307,shape=3, scale=1)

medianstandardizedmoments(x)

standardizedmoments(x)

#if the distribution is not Weibull, small biases might exist, but much smaller than others.
x<-rnorm(1*10^6)

SWAmoments(x)

medianstandardizedmoments(x)

imoments(x)

istandardizedmoments(x)

standardizedmoments(x)

#even for the kurtosis, the bias is relatively small, in fact, this can be further improved, but due to length limit, this is all at this stage.
