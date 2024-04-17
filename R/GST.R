#' GST function
#'
#' This function performs one or two-sided non-parametric global statistical test
#'
#' @param data a matrix of raw data. Rows are subjects and columns are endpoints.
#'  All columns must be coded so that larger value is a better outcome.
#'  The first m rows in data are data from treatment 1, and the rest rows are data from treatment 2
#' @param m number of subjects in treatment 1
#' @param test.side 1 (one-sided test) or 2 (two-sided test)
#'
#' @return `p.value`: p value of the global statistical test
#' @return `h`: adjusting factor to O'Brien's test
#' @return `theta.hat`: a vector of treatment effect on individual outcomes
#' @return `GTE`: global treatment effect of all outcomes
#' @return `sigmaGTE`: standard deviation of GTE estimate
#'
#' @import stats
#'
#' @export

GST<- function(data, m, test.side) {
  N<- dim(data)[1]
  n<- N-m
  k<- dim(data)[2]

  ind.y<- data*0
  ind.y[(m+1):N,]<- 1  # indicate observations of y.

  R.matrix<- apply(data, 2, rank) #rank by column
  R<- apply(R.matrix * ind.y, 2, sum)#rank sum in group y
  R.s<- apply(R.matrix, 1, sum)

  # compute theta.hat
  theta.hat<- 2*(R-n*(N+1)/2)/m/n    ### theta.hat=W in the paper

  A1<- array(NA, dim=c(m,k))
  for( i in 1:m)
  {
    y1<- data[c(i, (m+1):N),]
    Ry<- apply(y1, 2, rank)
    A1[i,]<- 2*Ry[1,]-2-n +n*theta.hat
  }

  B1<- array(NA, dim=c(n,k))
  for( i in (m+1):N)
  {
    x1<- data[c(i,1:m),]
    Rx<- apply(x1, 2, rank)
    B1[i-m,]<- 2*Rx[1,]-2-m -m*theta.hat
  }

  x1<- data[1:m,]
  Rx<- apply(x1,2,rank)
  A2<- 2*Rx-1-m
  y1<- data[(m+1):N,]
  Ry<- apply(y1,2,rank)
  B2<- 2*Ry-1-n

  h<- N^2/m/n*sum(t(A1)%*% A1+ t(B1)%*% B1)/
    sum(t(A1+A2)%*% (A1+A2)+t(B1+B2)%*% (B1+B2))

  # compute p value for rank-sum-type test
  if (test.side==1)
  {A<- t.test(R.s[(m+1):N],R.s[1:m],
              alternative="greater")
  p.value<-  1*(1-pt(abs(A$statistic/sqrt(h)),A$parameter))
  }
  if (test.side==2)
  {A<- t.test(R.s[(m+1):N],R.s[1:m],
              alternative="two.sided")
  p.value<-  2*(1-pt(abs(A$statistic/sqrt(h)),A$parameter))
  }
  p.value<- as.numeric(p.value)
  # compute standard deviation of GTE
  x<- GST.parameter(data, m)
  sigmaGTE<- x$sigmaGTE #=standard deviation of GTE

  list(p.value=p.value, h=h, theta.hat=theta.hat,
       GTE=mean(theta.hat), sigmaGTE=sigmaGTE)
}
