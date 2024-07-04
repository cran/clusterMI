#' Tune the number of iterations for variable selection using varselbest
#'
#' @param res.varselbest an output from the varselbest function
#' @param linewidth a numerical value setting the widths of lines
#' @param linetype what type of plot should be drawn
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param graph a boolean. If FALSE, no graphics are ploted. Default value is TRUE
#' @param title the main title
#' @importFrom ggplot2 aes ggplot geom_line labs
#' @importFrom reshape2 melt
#' @references Bar-Hen, A. and Audigier, V. An ensemble learning method for variable selection: application to high dimensional data and missing values. ArXiv e-prints <arXiv:1808.06952>
#' @keywords internal


chooseB.intern <- function(res.varselbest, linewidth = 1, linetype = "dotdash", xlab = "B", ylab = "Proportion", graph = TRUE, title = NULL){
  id<-proportion <- variable <-NULL;#for CRAN submission only
  gridB.intern <- seq(1,length(res.varselbest$res.detail))
  
  #extraction des variables retenues et selectionnees
  res.select<-lapply(res.varselbest$res.detail,"[[","res.select")
  
  #calcul de la proportion de selection pour les differentes simu
  matprop <- (apply(sapply(res.select,"[[","SousSelect"),
                  1,cumsum)/apply(sapply(res.select,"[[","garde"),1,cumsum))
  matprop[is.nan(matprop)] <- 0
  colnames(matprop) <- names(res.varselbest$res$proportion)
  matprop.plot <- as.data.frame(matprop)
  matprop.plot$id <- gridB.intern
  plot_data <- melt(matprop.plot, id.var = "id", 
                    value.name = "proportion")
  res.ggplot <- ggplot(plot_data, aes(x = id, y = proportion, group = variable, colour = variable)) + 
    geom_line(linetype = linetype, linewidth = linewidth) + 
    labs(x = xlab, y = ylab)+ labs(title = title)
  if(graph){
    print(res.ggplot)
  }
  res.out <- list(matprop = matprop, res.ggplot = res.ggplot)
  return(res.out)
}
