#' Attribute informed brain connectivity
#'
#' An MCMC algorithm for fitting the joint network and item response theory (JNIRT) model
#' to network data as well as item responses
#'
#' This command provides posterior inference for parameters in JNIRT models
#' assuming normal or binary data types
#'
#' "nrm": assuming the data is normally distributed, can be network or item responses
#'
#' "bin": assuming the data is binary, can be network or item responses
#'
#' @usage jnirt(X, Y, family_network, family_responses, K=0, rvar = FALSE ,
#' cvar = FALSE,  dcor = FALSE, model = UU,
#' intercept=TRUE, seed = 1, nscan =
#' 10000, burn = 500, odens = 25, plot=TRUE, print = TRUE, gof=TRUE,
#' prior=list())
#' @param X a list of V x V brain connectivity data.
#' @param Y a list of V x P attribute data.
#' @param W a matrix of N x Q covariates for the connectivity data.
#' @param H a matrix of N x Q1 covariates for the attribute data.
#' @param D integer: dimension of the multiplicative effects (can be zero) in the connectivity data.
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @param prior list: A list of hyperparameters for the prior distribution
#' @return
#' \item{VC}{posterior samples of the variance parameters}
#' \item{COV}{posterior samples of the covariance parameters}
#' \item{BETAPM}{posterior mean of the regression coefficient parameters for the connectivity data}
#' \item{GAMMAPM}{posterior mean of the regression coefficient parameters for the attribute data}
#' \item{THETAPM}{posterior samples of the latent person variable}
#' \item{APM}{posterior mean of connectivity intercepts} \item{BPM}{posterior
#' mean of attribute intercepts} \item{U}{posterior estimates of multiplicative
#' row effects u} \item{U_1}{the last iteration of the multiplicative
#' row effects u}
#' \item{UVPM}{posterior mean of UV}
#' \item{Theta_1}{the last iteration of the Theta estimate}
#' \item{X}{observed X}
#' \item{Y}{observed Y}
#' \item{UVC}{posterior samples of elements of the connectivity UU}
#' \item{TAC}{posterior samples of elements of the Theta}
#' \item{EFlPM}{posterior mean estimates of X}
#'  \item{ETPM}{posterior mean estiamtes of Y}
#' values of four goodness-of-fit statistics}
#' \item{input}{input values}
#' @author Selena Wang
#' @examples
#'
#' @export abc
#'




jbc_scaleAlpha_newU0_all<- function(X, Y,W, H, D=1,theoretical.str,
                    indices = NULL, indices_irt = NULL,
                    seed = 1, nscan = 10000, burn = 500, odens = 25,
                    print = TRUE, gof=TRUE, plot=TRUE,
                    prior=list())
{
  ## record model set up

  input<-list(D=D, theoretical.str,nscan=nscan, burn=burn, odens=odens, prior=prior, indices = indices, indices_irt = indices_irt)
  # set random seed
  set.seed(seed)

  K=1



  # starting Fl values
  Fl<-X
  Tm<-Y

  N<-length(X)
  V<-nrow(X[[1]])




  if(is.null(prior$Sutheta0)){ prior$Sutheta0<-diag(D+V) }
  if(is.null(prior$etautheta)){ prior$etautheta<-(D+V+2) }

  # starting intercept values
  a<-matrix(sapply(Fl, mean, na.rm=TRUE ))
  b<-matrix(colMeans(Tm,na.rm=TRUE))
  ET_b <- Tm- as.numeric(b)
  ET_b[is.na(ET_b)]=0


  # starting beta values
  if(!is.null(W)){ beta<-matrix(rep(0,ncol(W)), nrow = ncol(W),ncol=1) }else{beta<-NULL}
  if(!is.null(H)){ gamma<-matrix(rep(0,ncol(H)), nrow = ncol(H),ncol=1) } else{gamma<-NULL}



  ## Tm is behavior, s1 is behavior
  s1<-mean(sapply(1:length(X), function(x) mean((ET_b)^2, na.rm=TRUE)), na.rm=TRUE)
  s2<-mean(sapply(1:length(X), function(x) mean((Fl[[x]] - a[x,])^2, na.rm=TRUE)), na.rm=TRUE)




  # U
  tmp<-sapply(1:length(X), function(x) Fl[[x]] - a[x,], simplify = FALSE)

  for(i in 1:length(X)){tmp[[i]][is.na(tmp[[i]])]=0}




  Ul<-sapply(1:length(X), function(x) (svd(tmp[[x]])$u[,1:K,drop=FALSE]%*%diag(sqrt(svd(tmp[[x]])$d[1:K]),nrow=K)), simplify = FALSE)

  U=simplify2array(Ul)

  dimnames(U)<- list(c(seq(1,nrow(X[[1]]))),
                       c(1),
                       c(rownames(Y)))

  # Theta and Alpha

  Alpha<-matrix(0,ncol(Y),D)
  Theta<-matrix(0,nrow(Y),D)
  if(D>0)
  {
    sET<-svd(ET_b)
    Theta<-sET$u[,1:D,drop=FALSE]%*%diag(sqrt(sET$d[1:D]),nrow=D)
    Alpha<-sET$v[,1:D,drop=FALSE]%*%diag(sqrt(sET$d[1:D]),nrow=D)
  }






  # output items

  if(!is.null(W)){  BETA<- matrix(0,nrow = 0, ncol = ncol(W))}else{
    BETA<- matrix(0,nrow = 0, ncol = 0)
  }
  if(!is.null(H)){  GAMMA<- matrix(0,nrow = 0, ncol = ncol(H))}else{
    GAMMA<- matrix(0,nrow = 0, ncol = 0)
  }
  BETAPS <- beta * 0
  GAMMAPS <- gamma * 0

  ## behavior parameters
  bPS<- matrix(nrow = 0, ncol = ncol(Y))
  ALPHATHETAPS <- Theta %*% t(Alpha) * 0
  ALPHAPS <- Alpha*0
  THETAPS <- Theta*0
  UPS <-U*0
  XPM<-YPM<-EFlPM<-ETPM<-list()

  UVPS.l=sapply(1:(dim(U)[3]), function(i) U[,,i] %*% t(U[,,i]) * 0, simplify=FALSE )
  APS<-rep(0,length(X))
  names(APS)<- names(X)
  rownames(U)<-rownames(Theta)<-rownames(X[[1]])
  #XPS<-YPS<-list()

  #GOF<-list()
  #for(i in 1:length(X)){
  #gofXY<-c(gofstats_c(X[[i]]), gofstats_a(Y[[i]]))
  #GOF[[i]]<-matrix(gofXY,1,length(gofXY))
  #rownames(GOF[[i]])<-"obs"
  #colnames(GOF[[i]])<-names(gofXY)

  #YPS[[i]]<-matrix(0,nrow=V,ncol=P) ; dimnames(YPS[[i]])<-dimnames(Y[[i]])
  #XPS[[i]]<-matrix(0,nrow=V,ncol=V) ; dimnames(XPS[[i]])<-dimnames(X[[i]])

  #}



  if(is.null(indices)){
    indices<-matrix(sample(1:nrow(X[[1]]),min(round(nrow(X[[1]])/5),5)*2, replace = FALSE), nrow=2)

  }
  names_n<-NULL
  for(i in 1:ncol(indices)){names_n<-c(names_n,paste("UV",indices[1,i], indices[2,i],sep = ","))}

  if(is.null(indices_irt)){
    indices_irt<-matrix(c(sample(1:nrow(Y),min(round(nrow(X[[1]])/5),5), replace = FALSE),
                          sample(1:ncol(Y),min(round(nrow(X[[1]])/5),5), replace = TRUE)), nrow=2, byrow = TRUE)
  }

  names_i<-NULL
  for(i in 1:ncol(indices_irt)){names_i<-c(names_i,paste("ThetaAlpha",indices_irt[1,i], indices_irt[2,i],sep = ","))}
  TAC<-matrix(nrow=0,ncol=ncol(indices_irt))
  UVC<-matrix(nrow=0,ncol=ncol(indices))

  colnames(UVC) <- names_n
  colnames(TAC) <- names_i


  VC<-matrix(nrow=0,ncol=2+length(c(seq(1,(V)*(V+1)/2), seq(1,(D)*(D+1)/2))))

  COV<-matrix(nrow = 0, ncol = V*D)

  colnames(VC) <- c(paste("Su",seq(1,(V)*(V+1)/2),sep=""),
                    paste("Stheta",seq(1,(D)*(D+1)/2),sep=""), "ve_connectivity","ve_attributes")
  colnames(COV) <- paste("Sutheta",seq(1,(D)*V),sep="")

  ## prior for item parameters


  if(is.null(prior$mxi)){ prior$mxi<-c(rep(1,D),0) }
  if(is.null(prior$Sigmaxi)){ prior$Sigmaxi<-diag(D+1) }






  # MCMC
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))

  for (s in 1:(nscan + burn))
  {

    #check this
    # update Tm

    #if(is.null(H)){ET<-sapply(1:length(Y), function(x) (as.numeric(b[x,]) + Theta), simplify = FALSE)}
    #if(!is.null(H)){ET<-sapply(1:length(Y), function(x) (as.numeric(b[x,]) + as.numeric(H[x,] %*% gamma) + Theta), simplify = FALSE)}



    # update ET

    if(is.null(H)){ET <- matrix(rep(1,nrow(Y))) %*% t(b)+ Theta%*%t(Alpha)}
    if(!is.null(H)){ET <- matrix(rep(1,nrow(Y))) %*% t(b)+ Theta%*%t(Alpha) + H %*% gamma %*% matrix(rep(1,nrow(H)))^T}



    #if(family_responses=="nrm"){ H<-rH_nrm(H, EH,s1,Y) }
    #if(family_responses=="bin"){ H<-rH_bin(H, EH,Y) }



    # update Fl, U is an array, FL vs EFl, EFl is updated and Fl is not
    if(!is.null(W)){EFl<-sapply(1:length(X), function(x) (as.numeric(a[x,]) + as.numeric(W[x,] %*% beta) + U[,,x] %*% t(U[,,x])), simplify = FALSE)}
    if(is.null(W)){EFl<-sapply(1:length(X), function(x) (as.numeric(a[x,]) + U[,,x] %*% t(U[,,x])), simplify = FALSE)}


    # EFl=list()
    # for(i in 1:length(X)){
    #   if(!is.null(W)){EFl[[i]] <- as.numeric(a[i,]) + as.numeric(W[i,] %*% beta) + U %*% t(U)}else{
    #     EFl[[i]] <- as.numeric(a[i,]) + U %*% t(U)
    #   }
    #   #Fl[[i]] <- rFl_nrm(Fl[[i]], EFl[[i]],s2,X[[i]])
    #
    # }


    # update s2/s1 s1 is behavior
    s2<-rs2(Fl,offset = EFl, nu2=prior$nu2,s20=prior$s20)
    s1<-rs1_b(Tm,offset = ET, nu1=prior$nu1,s10=prior$s10)



    tmp<-rbeta_a_fc_per_b(Fl,W=W,s2=s2,U=U,ivA=prior$ivA,beta0=prior$beta0,S0=prior$S0)
    beta <- tmp$beta

    a<-tmp$a

    # update gamma,  b come back omitted for now
    #tmp1=rgamma_b_fc_per(ET,H=H,s1=s1,offset=Theta,ivB=prior$ivB,gamma0=prior$gamma0,S1=prior$S1)
    gamma <- NULL


    tmp <- rSu(data.matrix(cbind(t(U[,1,]),Theta)), Su0=prior$Sutheta0,etau=prior$etautheta)
    Su = matrix(tmp$Su[1:V,1:V], nrow=V, ncol=V)
    Stheta = matrix(tmp$Su[(V+1):(V+D),(V+1):(V+D) ], nrow=D, ncol=D)
    Sutheta =matrix(tmp$Su[(V+1):(V+D),1:V ], nrow = D, ncol = V)
    S=tmp$Su



    U<-rU_b_updated(Fl,U,Theta, Stheta, Sutheta, Su, s2, offset=sapply(1:length(EFl), function(x) EFl[[x]]-U[,,x]%*%t(U[,,x]), simplify = FALSE))

    U[,1,]=t(scale(t(U[,1,]), scale=FALSE))
    
    U[V,1,]=0

      # 
      # if(U[1,1,1]>0){
      #   U[,1,]=U[,1,]
      # }
      # 
      # if(U[1,1,1]<0){
      #   U[,1,]=U[,1,]*(-1)
      # }
      # 
      # 
    
    


    # if(s ==1){U_target<-U[,1,]}
    #
    # if(s>1){
    #   tmp <- Procrustes(U[,1,], U_target,
    #                     translate = FALSE,
    #                     dilate = FALSE,
    #                     sumsq = FALSE)
    #   U[,1,]<-tmp$X.new
    #
    #
    #
    # }


    ########
    Theta <-rTheta_b(Tm  , b, t(Alpha), t(Theta), t(U[,1,]), Stheta, Sutheta, Su, s1)
    #Theta <- scale(Theta)

    Theta <- Theta -   colMeans(Theta)



    #Theta <- t(mhalf(solve(cov(Theta))) %*% t(Theta))


    ####

    tmp.2 <-rXi(Tm, b, Alpha, Theta, prior$mxi, prior$Sigmaxi, s1)
    b <-tmp.2$beta
    Alpha_tmp <-tmp.2$Alpha

    tmp <- Procrustes(Alpha_tmp, theoretical.str,
                      translate = FALSE,
                      dilate = FALSE,
                      sumsq = FALSE)
    Alpha=Alpha_tmp %*% tmp$R
    Theta=Theta %*% tmp$R

    Alpha=matrix(1,nrow = 1,ncol=1)









    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn)
    {
      # save results

      VC<-rbind(VC, c( Su[upper.tri(Su, diag = T)],
                       Stheta[upper.tri(Stheta, diag = T)], s2, s1))
      COV<-rbind(COV, as.vector(Sutheta))
      ## covariates
      BETAPS<-BETAPS+beta
      GAMMAPS<-GAMMAPS+gamma


      UVPS.l <- sapply(1:(dim(U)[3]), function(i) UVPS.l[[i]] + U[,,i] %*% t(U[,,i]), simplify=FALSE)
      # update posterior sums of random effects
      #UVPS.l <- UVPS + U %*% t(U)

      tmp<-U[,,1] %*% t(U[,,1])
      tmp.1<-NULL
      for(i in 1:ncol(indices)){tmp.1 <- c(tmp.1 ,tmp[indices[1,i],indices[2,i]])}
      UVC<-rbind(UVC, c(tmp.1))
      ## behavior parameters
      THETAPS <- THETAPS + Theta
      ALPHAPS <- ALPHAPS + Alpha
      UPS<-UPS+U
      ALPHATHETAPS <- ALPHATHETAPS + Theta %*% t(Alpha)
      bPS<-rbind(bPS, as.vector(b))


      tmp<-Theta %*% t(Alpha)
      tmp.1<-NULL
      for(i in 1:ncol(indices_irt)){tmp.1 <- c(tmp.1 ,tmp[indices_irt[1,i],indices_irt[2,i]])}
      TAC<-rbind(TAC, c(tmp.1))

      ## intercept parameters
      APS <- APS + a
      #
      # Xs<-list()
      # Ys<-list()
      # for (i in 1:length(X)){
      #   Xs[[i]]<-simX_nrm(EFl[[i]],s2)
      #   Ys[[i]]<-simY_nrm(ET[[i]],s1)
      #   # update posterior sum
      #   XPS[[i]]<-XPS[[i]]+Xs[[i]]
      #   YPS[[i]]<-YPS[[i]]+Ys[[i]]
      #
      #   # save posterior predictive GOF stats
      #   #if(gof){ GOF[[i]]<-rbind(GOF[[i]],c(gofstats_c(Xs[[i]]),gofstats_a(Ys[[i]]))) }
      #
      # }
    }


  } # end MCMC

  # output


  # posterior means
  ### covariates
  GAMMAPM<-GAMMAPS/nrow(VC)
  BETAPM<-BETAPS/nrow(VC)

  ## intercept
  APM<-APS/nrow(VC)
  ## behavior parameters
  bPM<-colMeans(bPS)


  UVPM.l <- sapply(1:length(UVPS.l), function(i) UVPS.l[[i]]/nrow(VC), simplify=FALSE)


  ##
  ALPHAPM <-ALPHAPS/nrow(VC)
  THETAPM <-THETAPS/nrow(VC)
  UPM <-UPS/nrow(VC)

  # XPM<-sapply(XPS, function(x) x/nrow(VC), simplify = FALSE)
  # YPM<-sapply(YPS, function(x) x/nrow(VC), simplify = FALSE)


  if(!is.null(W)){EFlPM<-sapply(1:length(X), function(x) (APM[x] + as.numeric(W[x,] %*% BETAPM) + UVPM.l[[x]]), simplify = FALSE)}
  if(is.null(W)){EFlPM<-sapply(1:length(X), function(x) (APM[x] + UVPM.l[[x]]), simplify = FALSE)}




  if(!is.null(H)){ETPM<-matrix(rep(1,nrow(Y))) %*% t(bPM)+ THETAPM%*%t(ALPHAPM) + H %*% GAMMAPM %*% matrix(rep(1,nrow(H)))^T}
  if(is.null(H)){ETPM<-matrix(rep(1,nrow(Y))) %*% t(bPM)+ THETAPM%*%t(ALPHAPM)}


  # for(i in 1:length(X)){
  #   XPM[[i]]<-XPS[[i]]/nrow(VC)
  #   YPM[[i]]<-YPS[[i]]/nrow(VC)
  #
  #   if(!is.null(W)){ EFlPM[[i]]<-APM[i,] + as.numeric(W[i,] %*% BETAPM) + UVPM }else(
  #     EFlPM[[i]]<-APM[i,] + UVPM
  #   )
  #   if(!is.null(H)){ ETPM[[i]]<-BPM[i,] + as.numeric(H[i,] %*% GAMMAPM) + THETAPM }else{
  #     ETPM[[i]]<-BPM[i,] + THETAPM
  #   }
  #
  # }
  #

  names(APM)<-names(X)
  names(UVPM.l)<-rownames(THETAPM)<-rownames(X[[1]])
  rownames(ALPHAPM)<-colnames(Y)

  rownames(U)<-rownames(X[[1]])

  fit <- list(BETAPM=BETAPM,GAMMAPM=GAMMAPM, VC=VC,COV=COV, APM=APM,bPM=bPM,U=U,UVPM.l=UVPM.l,THETAPM = THETAPM,UPM=UPM, ALPHAPM = ALPHAPM,EFlPM=EFlPM,ETPM=ETPM,
              X=X,Y=Y, UVC=UVC, TAC=TAC, Theta=Theta, Alpha=Alpha, b=b, input=input, indices=indices, indices_irt=indices_irt)

  class(fit) <- "jbc"
  fit
}



