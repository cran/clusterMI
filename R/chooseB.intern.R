#' Tune the number of iterations for variable selection using varselbest
#'
#' @param res.varselbest an output from the varselbest function
#' @param gridB a grid for the number of iterations. By default, the grid is tuned to 1:B where B is the argument used in varselbest.
#' @param xlim the x limits (x1, x2) of the plot
#' @param plotvar index of variables for which a curve is ploted
#' @param cex a numerical vector giving the amount by which plotting characters and symbols should be scaled relative to the default
#' @param type what type of plot should be drawn
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param pch an integer
#' @param graph a boolean. If FALSE, no graphics are ploted. Default value is TRUE
#' @references Bar-Hen, A. and Audigier, V. An ensemble learning method for variable selection: application to high dimensional data and missing values. ArXiv e-prints <arXiv:1808.06952>
#' @keywords internal


chooseB.intern<-function(res.varselbest,gridB=NULL,xlim=NULL,plotvar=NULL,cex=.2,type="b",xlab="B",ylab="proportion",graph=TRUE, pch=16){
  if(is.null(gridB)){gridB.intern<-seq(1,length(res.varselbest$res.detail),ceiling(length(res.varselbest$res.detail)/100))}else{gridB.intern<-gridB}
  if(is.null(plotvar)){plotvar.intern<-seq(length(res.varselbest$res$proportion))}else{plotvar.intern<-plotvar}
  if(is.null(xlim)){xlim.intern<-gridB.intern}else{xlim.intern<-xlim}
  
  #extraction des variables retenues et s?lectionn?es
  res.select<-lapply(res.varselbest$res.detail,"[[","res.select")
  
  #calcul de la proportion de selection pour les diff?rentes simu
  matprop<-(apply(sapply(res.select,"[[","SousSelect"),
                  1,cumsum)/apply(sapply(res.select,"[[","garde"),1,cumsum))
  matprop[is.nan(matprop)]<-0
  
  if(graph){
    matplot(gridB.intern,matprop[gridB.intern,plotvar.intern],
            cex=cex,type=type,xlab=xlab,ylab=ylab, pch=pch)
  }
  
  return(matprop)
}
