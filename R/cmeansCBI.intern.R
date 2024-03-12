#' @importFrom e1071 cmeans

cmeansCBI.intern <- function (data, k = NULL, scaling = FALSE, m.cmeans=2) 
{
  if (!is.null(k)) 
    krange <- k
  if (!identical(scaling, FALSE)) 
    sdata <- scale(data, center = TRUE, scale = scaling)
  else sdata <- data
  c1 <- cmeans(x=sdata, centers=krange, m=m.cmeans)
  partition <- c1$cluster
  cl <- list()
  nc <- k
  for (i in 1:nc) cl[[i]] <- partition == i
  out <- list(result = c1, nc = k, clusterlist = cl, partition = partition, 
              clustermethod = "cmeans")
  out
}
