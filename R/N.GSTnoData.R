#' N.GSTnoData function
#'
#' This function calculates sample size when no primary data is available
#'
#' @param k number of endpoints
#' @param r0 n2/n1 randomization ratio
#' @param theta theta is a positive GTE where power is controlled
#' @param sigma2 upper bound of variances F1(X2) and F2(X1)
#' @param rho upper bound of absolute correlation between F1(X2) and F2(X1)
#' @param alpha type I error
#' @param power desired statistical power
#' @param test.side 1(one-sided test) or 2 (two-sided test)
#'
#' @return sample size needed from treatment 1
#' @return sample size needed from treatment 2
#' @return total sample size
#'
#' @import stats
#'
#' @export

N.GSTnoData<- function(k, r0, theta, sigma2, rho, alpha, power, test.side) {
  z.alpha<- qnorm(1-alpha/test.side)
  z.beta<- qnorm(power)
  N<- 4*sigma2/(theta^2) *(rho+1/k*(1-rho))*
    (1+r0)^2/r0*(z.alpha+z.beta)^2
  N<- round(N)
  n1<- round(N/(1+r0))
  n2<- round(N-n1)

  print(paste("sample size in treatment 1 =", n1))
  print(paste("sample size in treatment 2 =", n2))
  list(n1=n1, n2=n2, N=N)
}
