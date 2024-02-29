library(foreign)
library(dplyr)
library(gdata)
library(psych)
library(matrixcalc)

library(MASS)
library(expm)
library(boot)
library(Matrix)

library(CovTools)
library(ggplot2)
library(plyr)
library(lvm4net)
library(pROC)

library(sbm)

sourceEntireFolder <- function(folderName, verbose=FALSE, showWarnings=TRUE) {
  files <- list.files(folderName, full.names=TRUE)
  
  # Grab only R files
  files <- files[ grepl("\\.[rR]$", files) ]
  
  if (!length(files) && showWarnings)
    warning("No R files in ", folderName)
  
  for (f in files) {
    if (verbose)
      cat("sourcing: ", f, "\n")
    ## TODO:  add caught whether error or not and return that
    try(source(f, local=FALSE, echo=FALSE), silent=!verbose)
  }
  return(invisible(NULL))
}

## change it to code folder
sourceEntireFolder("/Users/selena/Desktop/github/LatentSNA/code", verbose=FALSE, showWarnings=TRUE)


#### generate data sample size
N<-50

## number of brain regions
V<-20
## number of behavior variables
P<-1




W<-NULL
H<-NULL

a_t<-matrix(0, nrow = N, ncol = 1)

D=1
## start to create covariance matrix
S <- diag(1,V+D)



n_signa=V


id_siga=c(6,3,12,17)

for (each in id_siga){
  
  S[each,id_siga[!id_siga %in% each]]=.9
}

S[(V+1),id_siga]=.9

S[id_siga,(V+1)]=.9




Su = matrix(S[1:V,1:V], nrow=V, ncol=V)
Stheta = matrix(S[(V+1):(V+D),(V+1):(V+D) ], nrow=D, ncol=D)
Sutheta =matrix(S[(V+1):(V+D),1:V], nrow = D, ncol = V)


UTheta <- mvrnorm(n = N, mu=rep(0,(V+D)), Sigma=S, empirical = FALSE)

U.array=array(NA, dim = c(V,1,N))

U.array[,1,]=t(UTheta[,1:(V)])
Theta_t <- data.matrix(UTheta[,(V+1):(V+D)])



## covariate coefficients

beta_t=NULL
gamma_t=NULL

#connectivity variance
#attribute variance
s1_t=.5





b_t=matrix(0, nrow = P, ncol = 1)

X<-list()
for(i in 1:N){
  
  errormatrix=matrix(0, nrow = V, ncol = V)
  errormatrix[upper.tri(errormatrix,diag = FALSE)]=rnorm(V*(V-1)/2, sd=sqrt(s2_t))
  errormatrix=t(errormatrix)+errormatrix
  diag(errormatrix)=rnorm(V, sd=sqrt(s2_t))
  
  X[[i]]<-as.numeric(a_t[i,])  + U.array[,,i]%*% t(U.array[,,i]) +errormatrix
  diag(X[[i]])=NA
  
}

Y= matrix(rep(1,N)) %*% t(b_t)+ Theta_t%*%t(Alpha_t)+matrix(rnorm(N*P, sd=sqrt(s1_t)),N,P)

### data stopped here




## fit model
model1=latentSNA(X, Y,W=NULL, H=NULL,
                   indices = NULL, indices_irt = NULL,
                   seed = 1, nscan = 1, burn = 1, odens = 1,
                   print = TRUE, gof=TRUE, plot=TRUE,
                   prior=list())





