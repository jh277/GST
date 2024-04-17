#' N.GSTwPilotData function
#'
#' @description
#' This function calculates sample size when some pilot data are available.
#'
#' @details
#' It returns the sample size needed in each treatment group based on
#' target treatment effect size, significance level, desired power for a specified
#' alternative treatment effect size, correlation among multiple outcomes, and finally,
#' the randomization ratio.
#'
#' @param pilot.data a (m1+m2) x k matrix of pilot data with its first m1 rows
#'  from treatment 1, and the rest m2 rows from treatment 2. k is the number of endpoints
#' @param m1 number of patients in treatment 1 in the pilot data
#' @param r0 n2/n1, randomization ratio in the study
#' @param test.side 1 (one-sided test) or 2 (two-sided test)
#' @param alpha type I error to be controlled
#' @param power statistical power to be controlled
#' @param GTE global treatment effect where power is controlled at.
#'  If GTE=NA (missing), the power will be controlled at the unbiased estimate
#'  of GTE (=mean of theta.pilot) from the pilot data.
#'
#' @return `n1`: sample size needed in treatment 1
#' @return `n2`: sample size needed in treatment 2
#' @return `theta.pilot`: unbiased estimate of GTE from the pilot data
#'
#' @import stats
#'
#' @export

N.GSTwPilotData<- function(pilot.data, m1, r0, test.side, alpha, power, GTE) {
  z.alpha<- qnorm(1-alpha/test.side)
  z.beta<- qnorm(power)
  x<- GST.parameter(pilot.data, m1)
  theta.pilot<- x$theta
  if (is.na(GTE)==T){GTE<- mean(theta.pilot)}
  m2<- nrow(pilot.data)-m1
  k<- ncol(pilot.data)
  d1.m1m2<- r0*(m1-1)*(m2-1)  #=d1/m1/m2
  d2<- (z.alpha+z.beta)^2/k^2/GTE^2
  # compute gamma1.m1m2=gamma1/m1/m2 and gamma2.m1m2=gamma2/m1/m2
  gamma1.m1m2<-(4*m1 + 4*m1^2*r0)*sum(x$S1.n1n2) +
    (4*m2^2 + 4*m2*r0)*sum(x$S2.n1n2) -
    (m1*m2^2 + m1^2*m2*r0)*sum(x$S3.n1n2)
  gamma2.m1m2<- 4*m1^2*sum(x$S1.n1n2) +
    4*m2^2*sum(x$S2.n1n2) - m1^2*m2^2*sum(x$S3.n1n2)
  # The required total sample size N is:
  N<- ceiling((1+r0)*d2*gamma1.m1m2/2/d1.m1m2*(1+sqrt(1-4*d1.m1m2*gamma2.m1m2/d2/gamma1.m1m2^2)))
  n1<- round(N/(1+r0))
  n2<- round(N-n1)

  print(paste("sample size in treatment 1 =", n1))
  print(paste("sample size in treatment 2 =", n2))
  print(paste("unbiased estimate of GTE =", mean(theta.pilot)))

  list(n1=n1, n2=n2, theta.pilot=theta.pilot)
}
