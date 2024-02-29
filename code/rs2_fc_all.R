#' Gibbs update for dyadic variance
#'
#' Gibbs update for dyadic variance
#'
#'
#' @usage rs2_fc(Z, rho,offset=0,nu0=NULL,s20=NULL)
#' @param Z n X n normal relational matrix
#' @param rho current value of rho
#' @param offset matrix of the same dimension as Z. It is assumed that Z-offset
#' is equal to dyadic noise, so the offset should contain any additive and
#' multiplicative effects (such as \code{Xbeta(X,beta+ U\%*\%t(V) +
#' outer(a,b,"+")  } )
#' @param nu0 prior degrees of freedom
#' @param s20 prior estimate of s2
#' @return a new value of s2
#' @author Peter Hoff
#' @export rs2_fc
rs2_fc_all <-
function(Fl,offset=0,nu2=NULL,s20=NULL)
{

    E=sapply(1:length(Fl), function(x) Fl[[x]]-offset[[x]], simplify = FALSE)
    length_all=sum(sapply(E, function(x) get_s2(x)$tmp.l, simplify = TRUE))
    sum_all=sum(sapply(E, function(x) get_s2(x)$tmp.s, simplify = TRUE))


    if(is.null(nu2)){ nu2<-1 }
    if(is.null(s20)){ s20<-1 }

    1/rgamma(1, (length_all+nu2)/2 , (sum_all+nu2*s20)/2 )
}
