#' Diagnostic plot for the number of iterations used in sequential imputation methods
#' 
#' @description The \code{choosemaxit} function plots the within and between variance for each variable (specified in \code{plotvars}) against the iteration number for each of the replications (specified in \code{plotm}).
#' 
#' @param output an outpout from the imputedata function
#' @param plotvars index of variables for which a curve is plotted
#' @param plotm a vector indicating which imputed datasets must be plotted
#' @param size size of points 
#' @param linewidth a numerical value setting the widths of lines
#' @param linetype what type of plot should be drawn
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param title the main title
#' @param nvar_by_row the number of variables that are plotted per window. Default value is 5.
#' @returns No return value
#' @export
#' @importFrom ggplot2 aes ggplot geom_line labs geom_point facet_grid
#' @importFrom reshape2 melt
#' @importFrom gridExtra arrangeGrob
#' @examples
#' data(wine)
#' set.seed(123456)
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- prodna(wine.na)
#' nb.clust <- 3 # number of clusters
#' m <- 3 # number of imputed data sets
#' maxit <- 50 # number of iterations for FCS imputation
#' \donttest{
#' res.imp <- imputedata(data.na = wine.na, method = "FCS-homo",
#'                       nb.clust = nb.clust, m = m, maxit = maxit)
#' choosemaxit(res.imp)
#' }


choosemaxit<-function(output,plotvars=NULL,plotm=1:5,size=.5,
                      linewidth = 1, linetype = "dotdash",xlab="iterations",ylab="var",
                      title= "Within and between variance plots",
                      nvar_by_row=5){
  m<-iter<-value<-id<-proportion <- variable <-NULL;#for CRAN submission only
  if(output$call$method%in%c("JM-DP","JM-GL")){stop("choosemaxit is dedicated to FCS imputation methods (FCS-homo or FCS-hetero)")}
  if(is.null(plotvars)){indexvar<-seq(dim(output$res.conv)[3])}else{indexvar<-plotvars}
  col<-plotm[which(plotm%in%seq(dim(output$res.conv)[1]))]
  res.conv.df<-reshape2::melt(output$res.conv)
  colnames(res.conv.df)<-c("m","iter","var","variance","value")
  res.conv.df <- res.conv.df[res.conv.df$m%in%plotm,,drop=FALSE]
  res.conv.df$m<-as.factor(res.conv.df$m)
  res.conv.df <-  droplevels(res.conv.df[res.conv.df$var%in%(levels(res.conv.df$var)[indexvar]),])
  levels(res.conv.df$variance) <- c("between","within")
  continue <- TRUE
  myplot <- list(); comp <- 1
  while(continue){
    myplot[[comp]]<-ggplot(res.conv.df[res.conv.df$var%in%unique(res.conv.df$var)[seq(from = (1+nvar_by_row*(comp-1)), to = nvar_by_row*(comp))],]) +
      aes(x = iter, y = value, colour = m) +
      geom_point(size=size) +
      geom_line(linewidth = linewidth, linetype = linetype, aes(x = iter, y = value, colour = m))+
      facet_grid(var~variance, scales = "free_y") + labs(x = xlab, y = ylab, title = title)
    comp <- comp+1
    continue <- (1+nvar_by_row*(comp-1))<length(unique(res.conv.df$var))
  }
  res.plot <- marrangeGrob(myplot,ncol = 1, nrow=1)
  print(res.plot)
  return()
}
