#' @keywords internal
randindex<-function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  RI <- (a + d)/(a + b + c + d)
  return(RI)
}