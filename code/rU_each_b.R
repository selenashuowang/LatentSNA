
rU_each_b<-function(u, n, Es,U,Theta, invUTheta, invThetaU, ivU, to, K,N){
  
  
  tmp.U=sapply(Es, function(x) apply(U[,,n]*x[u,],2,sum,na.rm=TRUE)-  U[u,,n]*x[u,u] , simplify = FALSE)
  #tmp.U1=Reduce('+', tmp.U)
  
  l<- tmp.U[[n]]*(to) - .5*t(invUTheta) %*% Theta[n,]- .5*invThetaU %*% Theta[n,]
  iQ<- solve( ivU +    (to)^2*( crossprod(U[,,n]) - U[u,,n]%*%t(U[u,,n]) ) ) 
  return(iQ%*%l + t(chol(iQ))%*%rnorm(K) )
}


