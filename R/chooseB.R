#' Diagnostic plot for the number of iterations used in the varselbest function
#'
#' @description 
#' \code{chooseB} plots the proportion of times an explanatory variable is selected according to the number of iterations (B).
#' @details
#' \code{varselbest} performs variable selection on random subsets of variables and, then, combines them to recover which explanatory variables are related to the response, following Bar-Hen and Audigier (2022) <doi:10.1080/00949655.2022.2070621>.
#' More precisely, the outline of the algorithm are as follows: let consider a random subset of \code{sizeblock} among p variables.
#' Then,  any selection variables scheme can be applied.
#' By resampling \code{B} times, a sample of size \code{sizeblock} among the p variables, we may count how many times a variable is considered as significantly related to the response and how many times it is not.
#' The number of iterations \code{B} should be large so that the proportion of times a variable is selected becomes stable. \code{chooseB} plots the values of proportion according to the number of iterations.
#' @param res.varselbest an output from the varselbest function
#' @param gridB a grid for the number of iterations. By default, the grid is tuned to \code{1:B} where \code{B} is the argument used in varselbest.
#' @param xlim the x limits (x1, x2) of the plot
#' @param plotvars index of variables for which a curve is plotted
#' @param cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default
#' @param type what type of plot should be drawn
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param graph a boolean. If FALSE, no graphics are plotted. Default value is TRUE
#' @param pch an integer
#' @return a list of matrices where each row corresponds to the vector of proportions (for all explanatory variables) obtained for a given value of B
#' @export
#' @references Bar-Hen, A. and Audigier, V., An ensemble learning method for variable selection: application to high dimensional data and missing values, Journal of Statistical Computation and Simulation, <doi:10.1080/00949655.2022.2070621>, 2022.
#' @seealso \code{\link{varselbest}}
#' @examples
#' data(wine)
#' 
#' require(parallel)
#' ref <- wine$cult
#' nb.clust <- 3
#' wine.na<-wine
#' wine.na$cult <- NULL
#' wine.na <- as.matrix(wine.na)
#' wine.na[sample(seq(length(wine.na)), size = ceiling(length(wine.na)/3))] <- NA
#' 
#' nnodes <- 2 # Number of CPU cores for parallel computing
#' B <- 80 # Number of iterations for variable selection
#' 
#' # variable selection
#' \donttest{
#' res.varsel <- varselbest(data.na = wine.na,
#'                         listvar = "alco",
#'                         B = B,
#'                         nnodes = nnodes,
#'                         nb.clust = nb.clust,
#'                         graph = FALSE)
#'# convergence
#'res.chooseB <- chooseB(res.varsel)
#'}

chooseB<-function(res.varselbest,gridB=NULL,xlim=NULL,plotvars=NULL,cex=.2,type="b",xlab="B",ylab="proportion",graph=TRUE, pch=16){
  res.out<-lapply(res.varselbest$res.varsel,
                  FUN = chooseB.intern,
                  gridB=gridB,
                  xlim=xlim,
                  plotvar=plotvars,
                  cex=cex,
                  type=type,
                  xlab=xlab,
                  ylab=ylab,
                  graph=graph,
                  pch=pch)
  return(res.out)
}
