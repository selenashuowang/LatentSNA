
rU_each_b_updated_abcd<-function(u, n, tmp.U.new.1,U,Theta, Stheta, Sutheta, Su, to, K,N){



  invUTheta <- (-1)*solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta)) %*% Sutheta %*% solve(Su)

  invThetaU <- (-1)* solve( Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta) %*% t(Sutheta) %*% solve(Stheta)

  ivU <- solve(Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta)
  
  
  tmp.new=matrix(c(U[-u,1,n],Theta[n,]))


  #tmp.U=sapply(Es, function(x) sum(U[,,n]*x[u,],na.rm=TRUE)-  U[u,,n]*x[u,u] , simplify = FALSE)
  #tmp.U1=Reduce('+', tmp.U)

  l<- tmp.U.new.1*(to) - .5*t(invUTheta) %*% tmp.new- .5*invThetaU %*% tmp.new
  iQ<- solve( ivU +    (to)^2*( crossprod(U[,,n]) - U[u,,n]%*%t(U[u,,n]) ) )
  return(iQ%*%l + t(chol(iQ))%*%rnorm(K) )
}


