#' Introduce missing values using a missing completely at random mechanism
#'
#' \code{prodna} generates an incomplete dataset by removing a proportion of observed values
#'
#' @return an incomplete data frame
#' 
#' @param X a data frame (or matrix).
#' @param pct proportion of missing values. By default pct = 0.3.
#' @export
#' @examples
#' n <- 1000
#' p <- 5
#' X <- matrix(runif(n*p), nrow = n, ncol = p)
#' summary(X)
#' X.na <- prodna(X)
#' colMeans(is.na(X.na))

prodna<-function(X, pct = 0.3){
  Xmat<-as.matrix(X)
  ismissing <- sample(seq.int(length(Xmat)), size = ceiling(length(Xmat)*pct))
  Xmat[ismissing] <- NA
  res.out <- as.data.frame(Xmat)
  return(res.out)
}
