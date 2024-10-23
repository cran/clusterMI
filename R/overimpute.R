#' Overimputation diagnostic plot
#' 
#' @description
#' \code{overimpute} assesses the fit of the predictive distribution after performing multiple imputation with the \code{\link{imputedata}} function
#' 
#' 
#' @details
#' This function imputes each observed value from each conditional imputation model obtained from the imputedata function. The comparison between the "overimputed" values and the observed values is made by building a confidence interval for each observed value using the quantiles of the overimputed values (see Blackwell et al. (2015) <doi:10.1177/0049124115585360>). Note that confidence intervals built with quantiles require a large number of imputations. If the model fits well the data, then the 90\% confidence interval should contain the observed value in 90\% of the cases. The function overimpute takes as an input an output of the \code{\link{imputedata}} function (\code{res.imputedata} argument), the indices of the incomplete continuous variables that are plotted (\code{plotvars}), the indices of individuals (can be useful for time consuming imputation methods), the number of CPU cores for parallel computation, and the path for exporting print message generated during the parallel process (\code{path.outfile}).
#' 
#' @return A list of two matrices
#'   \item{res.plot}{7-columns matrix that contains (1) the variable which is overimputed, (2) the observed value of the observation, (3) the mean of the overimputations, (4) the lower bound of the confidence interval of the overimputations, (5) the upper bound of the confidence interval of the overimputations, (6) the proportion of the other variables that were missing for that observation in the original data, and (7) the color for graphical representation}
#'   \item{res.values}{a matrix with overimputed values for each cell. The number of columns corresponds to the number of values generated (i.e. the number of imputed datasets)}
#' 
#' @param res.imputedata an output from the imputedata function
#' @param plotvars column index of the variables overimputed
#' @param plotinds row index of the individuals overimputed
#' @param nnodes an integer indicating the number of nodes for parallel calculation. Default value is 2
#' @param path.outfile a vector of strings indicating the path for redirection of print messages. Default value is NULL, meaning that silent imputation is performed. Otherwise, print messages are saved in the files path.outfile/output.txt. One file per node is generated.
#' @param alpha	alpha level for prediction intervals
#' @param mfrow a vector of the form c(nr, nc)
#' @param mar a numerical vector of the form c(bottom, left, top, right)
#' @export
#' @importFrom mix em.mix imp.mix da.mix prelim.mix rngseed
#' @importFrom mclust Mclust
#' @importFrom mice mice complete
#' @importFrom micemd mice.impute.2l.2stage.norm mice.impute.2l.2stage.pmm mice.impute.2l.glm.norm mice.impute.2l.jomo
#' @importFrom grDevices heat.colors
#' @importFrom graphics legend segments
#' @importFrom stats quantile
#' @references Blackwell, M., Honaker, J. and King. G. 2015. A Unified Approach to Measurement Error and Missing Data: Overview and Applications. Sociological Methods and Research, 1-39. <doi:10.1177/0049124115585360>
#' @examples 
#' data(wine)
#' 
#' require(parallel)
#' set.seed(123456)
#' ref <- wine$cult
#' nb.clust <- 3
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- prodna(wine.na)
#' \donttest{
#' nnodes <- 2 # Number of CPU cores used for parallel computation
#' 
#' # Multiple imputation using m = 100 (should be larger in practice)
#' 
#' res.imp.over <- imputedata(data.na = wine.na,
#'                            nb.clust = nb.clust,
#'                            m = 100)
#' # Overimputation
#' 
#' ## overimputed variable
#' plotvars <- "alco" 
#' 
#' ## selection of 20 complete individuals on variable "alco"
#' plotinds <- sample(which(!is.na(wine.na[, plotvars])),
#'                     size = 20)
#' ## overimputation                   
#' res.over <- overimpute(res.imp.over,
#'                        nnodes = nnodes,
#'                        plotvars = plotvars,
#'                        plotinds = plotinds,
#'                        )
#' }

overimpute <- function (res.imputedata, plotvars = NULL, plotinds = NULL, nnodes = 2, 
          path.outfile = NULL, alpha = 0.1, mfrow=NULL,mar=c(5, 4, 4, 2) - 1.9) 
{
  m.intern <- sum(!sapply(res.imputedata$res.imp,is.null))

  if (m.intern < 100) {
    warning("The number of imputed data sets is too low to build confidence intervals according to the quantiles method. You should run imputedata with m over than 100.")
  }
  if (is.null(plotvars)) {
    Var <- seq(ncol(res.imputedata$call$data.na))
  }
  else {
    Var <- plotvars
  }
  if (is.null(plotinds)) {
    Ind <- seq(nrow(res.imputedata$call$data.na))
  }
  else {
    Ind <- plotinds
  }

  is.plot <- sapply(as.data.frame(res.imputedata$call$data.na[, Var, drop = FALSE]), is.numeric)
  is.plot[apply(res.imputedata$call$data.na[, Var, drop = FALSE], 2, function(xx) {
    length(table(xx)) == 2
  })] <- FALSE

  is.plot[colSums(is.na(res.imputedata$call$data.na[, Var, drop = FALSE])) == 
            0] <- FALSE

  don.plot <- as.matrix(res.imputedata$call$data.na[Ind, Var[is.plot], drop = FALSE])
  res.over <- matrix(NA, nrow = sum(!is.na(don.plot)), ncol = m.intern)
  cl <- parallel::makeCluster(nnodes, type = "PSOCK")
  parallel::clusterExport(cl, list("is.plot", "res.imputedata", 
                         "path.outfile", "Ind","imputedata","myem.mix","em.mix", "myimp.mix","imp.mix", "da.mix", "prelim.mix", "rngseed", "Mclust",
                         "mice", "complete",
                         "readData", "createModel", "multipleImp",
                         "rdirichlet","mice.impute.2l.2stage.norm","mice.impute.2l.2stage.pmm","mice.impute.2l.glm.norm","mice.impute.2l.jomo"), envir = environment())
  if (!is.null(path.outfile)) {
    parallel::clusterEvalQ(cl, sink(paste0(path.outfile, "/output", 
                                 Sys.getpid(), ".txt")))
  }

  res.over <- t(parallel::parSapply(cl, which(!is.na(don.plot)), FUN = function(jj, 
                                                                      is.plot, res.imputedata, Ind,verbose) {
    if(verbose){cat("cell number : ", jj, "\n")}
    sapply(seq(m.intern), FUN = function(m, ii, is.plot, 
                                                      res.imputedata, Ind) {
      don.over <- res.imputedata$res.imp[[m]]
      tmp <- as.matrix(don.over[Ind, names(which(is.plot))])
      tmp[ii] <- NA
      don.over[Ind, names(which(is.plot))] <- tmp
      res.imp.tmp <- try(imputedata(don.over,
                                    m = 1,
                                    maxit = 1, 
                                    verbose=FALSE,
                                    method = res.imputedata$call$method,
                                    method.mice = res.imputedata$call$method.mice, 
                                    predictmat = res.imputedata$call$predictmat,
                                    nb.clust=res.imputedata$call$nb.clust,
                                    L = 2,
                                    Lstart = 1,
                                    bootstrap=res.imputedata$call$bootstrap))
      
      if (inherits(res.imp.tmp,"try-error")) {
        warning(res.imp.tmp)
        res <- NA
      }
      else {
        res <- as.matrix(res.imp.tmp$res.imp[[1]])[Ind, names(which(is.plot))][ii]
      }
      return(res)
    }, ii = jj, is.plot = is.plot, res.imputedata = res.imputedata, Ind = Ind)
  }, is.plot = is.plot, res.imputedata = res.imputedata, Ind = Ind, verbose = (!is.null(path.outfile))))
  parallel::stopCluster(cl)
  
  res.plot <- t(apply(res.over, 1, function(x, alpha) {
    xbar <- mean(x,na.rm=TRUE)
    temp <- quantile(x, probs = c(alpha/2, 1 - alpha/2), 
                     na.rm = TRUE)
    binf <- temp[[1]]
    bsup <- temp[[2]]
    return(c(xbar = xbar, binf = binf, bsup = bsup))
  }, alpha = alpha))
  pct <- rep(rowMeans(is.na(res.imputedata$call$data.na[Ind, ])), ncol(don.plot))
  col <- cut(pct, c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.1))
  levels(col) <- c("blue", "green", heat.colors(3)[c(3, 
                                                     2, 1)])
  col <- as.character(col)
  res.over.value <- res.over
  res.over <- cbind.data.frame(var = unlist(mapply(FUN = rep, 
                                                   names(is.plot[is.plot]),
                                                   each = apply(!is.na(don.plot), 2, sum))),
                               trueval = don.plot[!is.na(don.plot)], 
                               res.plot, pct = pct[!is.na(don.plot)],
                               col = col[!is.na(don.plot)])
  colnames(res.over)[1] <- "var"
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if(is.null(mfrow)){
    Mfrow<-c(ceiling(sqrt(length(is.plot[is.plot]))),
                              ceiling(length(is.plot[is.plot])/ceiling(sqrt(length(is.plot[is.plot])))))
  }else{
    Mfrow<-mfrow
  }
  pct <- unlist(by(res.over,INDICES = res.over$var,function(res.over){
    round(100*mean(((res.over[, "trueval"] <= res.over[, "bsup"]) & (res.over[, "trueval"] >= res.over[, "binf"]))), 
                                                                          2)},simplify = FALSE))
  pct.rep<-pct[res.over$var]
  res.over$pct <- pct.rep
  par(mfrow = Mfrow, mar = mar)
  by(res.over, INDICES = res.over$var, FUN = function(xx) {
    pct <- xx$pct[1]
    plot(x = xx[, "trueval"], y = xx[, "xbar"], col = as.character(xx[, 
                                                                      "col"]), xlab = "observed values", ylab = "imputed values", 
         main = paste(xx[1, "var"], " (cov =", pct,"%)"), 
         ylim = c(min(xx[, "binf"], na.rm = T), max(xx[, 
                                                       "bsup"], na.rm = T)))
    abline(0, 1)
    segments(x0 = xx[, "trueval"], x1 = xx[, "trueval"], 
             y0 = xx[, "binf"], y1 = xx[, "bsup"], col = as.character(xx[, 
                                                                         "col"]))
    legend("bottomright", legend = c("0-0.2", "0.2-0.4", 
                                     "0.4-0.6", "0.6-0.8", "0.8-1"), col = c("blue", 
                                                                             "green", heat.colors(3)[c(3, 2, 1)]), bty = "n", 
           lty = 1, horiz = TRUE, cex = 1, lwd = 0.4)
  })
  return(list(res.plot = res.over, res.values = res.over.value))
}

