

get_s2<-function(E){
  rho=0
  H<-mhalf( solve(matrix(c(1,rho,rho,1),2,2)) )
  EM<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] ) %*%H
  ED<-diag(E)/sqrt(1+rho)

  tmp.s=sum(EM^2)+sum(ED^2)
  tmp.l=length(EM)+length(ED)

  return(list("tmp.s"=tmp.s, "tmp.l"=tmp.l))
  }
