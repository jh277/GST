#' GST.parameter function
#'
#' This function calculate GST p.value and estimate GTE
#'
#' @param data N-by-k matrix whose first n1 rows are observations from group 1
#' @param n1 number of observations in group 1
#'
#' @return unbiased estimator of theta
#' @return estimated standard deviation of GTE
#' @return estimated S1.n1n2
#' @return estimated S2.n1n2
#' @return estimated S3.n1n2
#' @return unbiased estimator of theta sigma
#' @return estimated h1, adjusted rank-sum-type test
#' @return estimated h2, adjusted rank-sum-type test
#'
#' @import stats
#'
#' @export

GST.parameter<- function(data, n1) {
  # no missing values are allowed
  data<- as.matrix(data)
  if (is.na(sum(data))==T){stop("NAs are not allowed in data")}

  N<- dim(data)[1]
  n2<- N-n1
  k<- dim(data)[2]
  Sigma.hat<- matrix(NA,ncol=k,nrow=k)
  theta<- rep(NA,k)

  ind.2<- data*0
  ind.2[(n1+1):N,]<- 1  # indicate observations in group 2.

  # compute the total ranks for group 2
  R<- apply(data, 2, rank)
  R<- apply(R * ind.2, 2, sum) # group 2 rank sums

  # compute an unbiased estimate of theta
  theta.hat<- (2*R-n2*(N+1))/n1/n2
  GTE<- mean(theta.hat)
  A1<- array(NA, dim=c(n1,k))
  B1<- array(NA, dim=c(n2,k))

  # To avoid numerical overflow, we do not compute S1, S2, and S3
  # but compute their ratio with n1*n2

  # compute matrix S1.n1n2=S1/n1/n2:
  x2<- array(NA, dim=c(n2,k))
  for( i in (n1+1):N)
  {
    x1<- data[c(i,1:n1),]
    Rx<- apply(x1, 2, rank)
    x2[i-n1,]<- Rx[1,]
    B1[i-n1,]<- 2*Rx[1,]-2-n1 -n1*theta.hat
  }
  # S1.n1n2<- var(x2, unbiased=FALSE)/n1/n2 (Splus code)
  S1.n1n2<- var(x2)/n1/n2

  # compute matrix S2.n1n2=S2/n1/n2:
  x2<- array(NA, dim=c(n1,k))
  for( i in 1:n1)
  {
    x1<- data[c(i,(n1+1):N),]
    Ry<- apply(x1, 2, rank)
    x2[i,]<- Ry[1,]
    A1[i,]<- 2*Ry[1,]-2-n2 +n2*theta.hat
  }
  S2.n1n2<- var(x2)/n1/n2

  # compute matrix S3.n1n2=S3/n1/n2:
  x1<- data[rep(1:n1, times=n2),]
  x2<- data[rep(((n1+1):N), each=n1),]

  x3<- 1*(x1<x2)-1*(x1>x2)
  S3.n1n2<- var(x3)/n1/n2

  # compute an unbiased estimate of Sigma
  Sigma.hat<- N/(n1-1)/(n2-1)*(S1.n1n2*n1+S2.n1n2*n2-S3.n1n2/4*n1*n2)
  if(sum(Sigma.hat)<=0)
  {# the following Sigma.hat is only a consistent estimate of Sigma
    # Since matrix S3 is ignored:
    Sigma.hat<- N/(n1-1)/(n2-1)*(S1.n1n2*n1+S1.n1n2*n2)
  }

  # compute h1, h2 for adjusted rank-sum-type test
  x1<- data[1:n1,]
  Rx<- apply(x1,2,rank)
  A2<- 2*Rx-1-n1
  y1<- data[(n1+1):N,]
  Ry<- apply(y1,2,rank)
  B2<- 2*Ry-1-n2

  h1<- N^2/n1/n2*sum(t(A1)%*% A1+ t(B1)%*% B1)/
    sum(t(A1+A2)%*% (A1+A2)+ t(B1+B2)%*% (B1+B2))

  h2<- N^2*sum(t(A1)%*% A1+ t(B1)%*% B1)/
    sum(n2^2*t(A1+A2)%*% (A1+A2)+ n1^2*t(B1+B2)%*% (B1+B2))
  # standard deviation of GTE:
  sigmaGTE<- 2*sqrt(sum(Sigma.hat)/N)/k

  list(theta=theta.hat, sigmaGTE=sigmaGTE,
       S1.n1n2=S1.n1n2, S2.n1n2=S2.n1n2, S3.n1n2=S3.n1n2,
       Sigma=Sigma.hat, h1=h1,h2=h2)
}
