#' Simulate missing values in a normal item response model
#' 
#' Simulates missing values under a item response model
#' 
#' 
#' @usage rH_nrm(H, EH,s1, Y)
#' @param H a square matrix, the current value of H
#' @param EH expected value of H
#' @param s1 item response variance
#' @param Y item response matrix
#' @return an item response matrix, equal to  at non-missing values
#' @author Selena Wang
#' @export rH_nrm
rH_nrm<-function(H,EH,s1,Y)
{
  HS<-simY_nrm(EH,s1)
  H[is.na(Y)]<-HS[is.na(Y)]  # this isn't quite right if there is asymmetric missingness. 
  H
}

