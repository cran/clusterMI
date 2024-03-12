#' @importFrom stats model.matrix
#' @keywords internal

class2ind <- function (x, drop2nd = FALSE) 
{
  if (!is.factor(x)) 
    stop("'x' should be a factor")
  y <- model.matrix(~x - 1)
  colnames(y) <- gsub("^x", "", colnames(y))
  attributes(y)$assign <- NULL
  attributes(y)$contrasts <- NULL
  if (length(levels(x)) == 2 & drop2nd) {
    y <- y[, 1]
  }
  y
}