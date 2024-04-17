#' theta function
#'
#' This function calculates global treatment effect
#'
#' @param x1 observations from group 1
#' @param x2 observations from group 2
#'
#' @return `theta1`: unbiased estimator of theta, which is the treatment effect on individual outcomes
#'
#' @import stats
#'
#' @export


# theta:
theta<- function(x1,x2)
{
  n1<- length(x1)
  n2<- length(x2)
  R2<- sum(rank(c(x1,x2))[(n1+1):(n1+n2)]) #=rank of x2
  theta1<- (2*R2-n2*(n1+n2+1))/n1/n2
  theta1
}
