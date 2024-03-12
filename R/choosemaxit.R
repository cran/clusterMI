#' Diagnostic plot for the number of iterations used in sequential imputation methods
#' 
#' @description The \code{choosemaxit} function plots the within and between variance for each variable (specified in \code{plotvars}) against the iteration number for each of the replications (specified in \code{plotm}).
#' 
#' @param output an outpout from the imputedata function
#' @param plotvars index of variables for which a curve is plotted
#' @param plotm a vector indicating which imputed datasets must be plotted
#' @param cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default
#' @param pch a vector of plotting characters or symbols
#' @param type what type of plot should be drawn
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param layout a vector of the form c(nrow, ncol)
#' @returns No return value
#' @export
#' @importFrom grDevices rgb col2rgb
#' @importFrom reshape2 melt
#' @importFrom lattice xyplot
#' @examples
#' data(wine)
#' set.seed(123456)
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- as.matrix(wine.na)
#' wine.na[sample(seq(length(wine.na)), size = ceiling(length(wine.na)/3))] <- NA
#' nb.clust <- 3 # number of clusters
#' m <- 3 # number of imputed data sets
#' maxit <- 50 # number of iterations for FCS imputation
#' \donttest{
#' res.imp <- imputedata(data.na = wine.na, method = "FCS-homo",
#'                       nb.clust = nb.clust, m = m, maxit = maxit)
#' choosemaxit(res.imp, cex = .7)
#' }


choosemaxit<-function(output,plotvars=NULL,plotm=1:5,cex=.3,
                      pch=16,type="b",xlab="iterations",ylab="var",
                      layout=NULL){
  m<-NULL#for CRAN submission only
  if(output$call$method%in%c("JM-DP","JM-GL")){stop("choosemaxit is dedicated to FCS imputation methods (FCS-homo or FCS-hetero)")}
  if(is.null(plotvars)){indexvar<-seq(dim(output$res.conv)[3])}else{indexvar<-plotvars}
  p <- length(indexvar)
  if(is.null(layout)){
  Mfrow <- c(min(p, 4), 1 + (p - 1)%/%4)
  }else{
    Mfrow <- layout
  }
  col<-plotm[which(plotm%in%seq(dim(output$res.conv)[1]))]
  res.conv.df<-reshape2::melt(output$res.conv)
  colnames(res.conv.df)<-c("m","iter","var","variance","value")
  res.xyplot<-lattice::xyplot(value ~ iter |variance+var,
                groups= m,
                data = res.conv.df,
                superpose=T,
                col=sapply(sapply(col,FUN=function(xx){t(col2rgb(xx))/255},simplify = FALSE),rgb,alpha=.5),
                type=c("b"),
                main="Convergence diagnostic",
                scales=list(y=list(relation="free")),
                cex=.1,lty=2,layout=Mfrow)
  print(res.xyplot)
  return()
}
