#' GST.simu function
#'
#' This function performs simulation using GST
#'
#' @param data N-by-k matrix whose first n1 rows are observations from group 1
#' @param m number of observations in group 1
#'
#' @return `p.GST`: p value of the global statistical test
#' @return `GTE`: global treatment effect of all outcomes
#' @return `p.univariate`: p value under univariate assumption
#'
#' @import stats
#'
#' @export

GST.simu<- function(data, m)
{
  N<- dim(data)[1] #=m+n
  k<- dim(data)[2]
  simu<- dim(data)[3]
  p.GST<- rep(NA, simu)
  p.univariate<- matrix(NA, nrow=simu, ncol=k)
  GTE<- rep(NA, simu)
  for (i in 1:simu)
  {
    x<- GST(data[,,i], m, 1) #one-sided test
    p.GST[i]<- x$p.value
    GTE[i]<- x$GTE
    for(k in 1:k)
    {
      x<- t.test(data[1:m, k, i], data[(m+1):N, k,i])
      p.univariate[i,k]<- x$p.value
    }

  }
  list(p.GST=p.GST, GTE=GTE, p.univariate=p.univariate)
}
