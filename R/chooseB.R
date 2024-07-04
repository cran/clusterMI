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
#' @param plotvar index of variables for which a curve is ploted
#' @param linewidth a numerical value setting the widths of lines
#' @param linetype what type of plot should be drawn
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param nrow argument of gtable. Default value is 2.
#' @param ncol argument of gtable. Default value is 2.
#' @param graph a boolean. If FALSE, no graphics are ploted. Default value is TRUE
#' @importFrom gridExtra marrangeGrob
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
#' wine.na <- prodna(wine.na)
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
#' res.chooseB <- chooseB(res.varsel)
#'}

chooseB<-function(res.varselbest,plotvar = NULL, linewidth = 1, linetype = "dotdash", xlab = "B", ylab = "Proportion", nrow = 2, ncol = 2, graph = TRUE){
  
  if(!is.null(plotvar)){
    res.varsel.intern<-res.varselbest$res.varsel[plotvar]
  }else if (is.null(plotvar)){
    res.varsel.intern<-res.varselbest$res.varsel
  }
  res.out<-mapply(res.varselbest = res.varsel.intern,
                  title = names(res.varsel.intern),
                  FUN = chooseB.intern,
                  MoreArgs = list(
                    linewidth=linewidth,
                    linetype =linetype,
                    xlab=xlab,
                    ylab=ylab,
                    graph=FALSE),SIMPLIFY = FALSE)
 
  if(graph){
    res.plot <- lapply(res.out, "[[","res.ggplot")
    res.plot <- marrangeGrob(res.plot,ncol = ncol, nrow=nrow)
    print(res.plot)
  }
  
  res.out <- lapply(res.out,"[[","matprop")
  return(res.out)
}
