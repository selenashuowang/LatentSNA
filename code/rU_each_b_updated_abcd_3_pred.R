
rU_each_b_updated_abcd_3_pred<-function(u,utheta.tmp,U,Theta, Stheta, Sutheta, Su, to, K,N){



  invUTheta <- (-1)*solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta)) %*% Sutheta %*% solve(Su)

  invThetaU <- (-1)* solve( Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta) %*% t(Sutheta) %*% solve(Stheta)

  ivU <- solve(Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta)
  
  

  #tmp.U=sapply(Es, function(x) sum(U[,,n]*x[u,],na.rm=TRUE)-  U[u,,n]*x[u,u] , simplify = FALSE)
  #tmp.U1=Reduce('+', tmp.U)

  

  l= (-0.5)*matrix(t(invUTheta) %*% utheta.tmp)- .5*matrix(invThetaU %*% utheta.tmp)
  iq_l=diag(1/((to)^2*(colSums(U[,1,]^2)-(U[u,1,]^2)) +rep(ivU,N) ))
  return(iq_l%*%l + t(chol(iq_l))%*%rnorm(N)) 
    


}


