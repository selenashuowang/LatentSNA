
rU_each_b_updated<-function(u, n, Es,U,Theta, Stheta, Sutheta, Su, to, K,N){

  Sutheta_small= Sutheta[,u]
  Su_small= Su[u,u]

  invUTheta <- (-1)*solve(Stheta - Sutheta_small %*% solve(Su_small) %*% t(Sutheta_small)) %*% Sutheta_small %*% solve(Su_small)

  invThetaU <- (-1)* solve( Su_small - t(Sutheta_small) %*% solve(Stheta) %*% Sutheta_small) %*% t(Sutheta_small) %*% solve(Stheta)


  ivU <- solve(Su_small - t(Sutheta_small) %*% solve(Stheta) %*% Sutheta_small)


  tmp.U=sapply(Es, function(x) sum(U[,,n]*x[u,],na.rm=TRUE)-  U[u,,n]*x[u,u] , simplify = FALSE)
  #tmp.U1=Reduce('+', tmp.U)

  l<- tmp.U[[n]]*(to) - .5*t(invUTheta) %*% Theta[n,]- .5*invThetaU %*% Theta[n,]
  iQ<- solve( ivU +    (to)^2*( crossprod(U[,,n]) - U[u,,n]%*%t(U[u,,n]) ) )
  return(iQ%*%l + t(chol(iQ))%*%rnorm(K) )
}


