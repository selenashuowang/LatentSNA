#' Gibbs sampling of U
#'
#' A Gibbs sampler for updating U.
#'
#' @usage rU(Fl,U,Theta, Stheta, Sutheta, Su, s2=1, offset=offset)
#' @param Fl a list of V X V normal relational matrix
#' @param U V X K matrix containing current value of U
#' @param Theta D X V current value of Theta
#' @param Stheta D X D covariance of Theta
#' @param Sutheta D X K covariance between U and Theta
#' @param Su K X K matrix containing covariance of U
#' @param s2 dyadic variance
#' @param offset a list of the same dimension as Fl. It is assumed that
#' Fl-offset follows a SRRM, so the offset should contain any multiplicative
#' effects (such as \code{U\%*\% t(U) } )
#' @return \item{U}{a new value of U}
#' @author Selena Wang
#'
#'
#' @export rU
rU_b_updated_abcd_2<-function(Fl,EFl, U,Theta, Stheta, Sutheta, Su, S, s2)
{
  K<-dim(U)[2] ; V<-dim(U)[1] ; N <-dim(U)[3]


  to<-as.numeric(sqrt(solve(s2)))

  Es=sapply(1:length(Fl), function(x) (Fl[[x]]-EFl[[x]]+U[,,x]%*%t(U[,,x]))*to, simplify = FALSE)
  Es=lapply(Es, function(x) { if(is.na(diag(x)[1])){diag(x) <- 0}; x})

  
  tmp.fz.new=sapply(1:V, function(u) sapply(1:N, function(n) sum(U[,,n]*Es[[n]][u,],na.rm=TRUE)-  U[u,,n]*Es[[n]][u,u] , simplify = TRUE),simplify = TRUE)
  utheta.tmp.list=sapply(1:V, function(u) sapply(1:N, function(n) matrix(c(U[-u,1,n],Theta[n,])) , simplify = TRUE),simplify = FALSE)
  
  #tmp.new=matrix(c(U[-u,1,n],Theta[n,]))
  return(sapply(1:V, function(u) rU_each_b_updated_abcd_2(u,  tmp.fz.new[,u],utheta.tmp.list[[u]],U,Theta,  Stheta=S[-u,-u], Sutheta=S[-u,u], Su=S[u,u], to, K,N), simplify = TRUE))



  # start=Sys.time()
  # for (n in 1:N){
  #   U[,,n]=matrix(sapply(1:V, function(x) rU_each_b_updated(u=x, n, Es,U,Theta,  Stheta, Sutheta, Su, to, K,N), simplify = TRUE))
  # }
  # 
  # print(Sys.time()-start)
  # 

}







