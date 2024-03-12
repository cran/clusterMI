#' MclustBootstrap with nboot = 1 and the same output as Mclust
#'
#' @param don matrix, or data frame
#' @param G An integer vector specifying the numbers of mixture components
#' @param modelNames A vector of character strings indicating the models to be fitted 
#' @param verbose A logical controlling if a text progress bar is displayed
#' @keywords internal
#' @importFrom mclust Mclust hc

mclustboot.intern<-function(don, G = NULL, modelNames = NULL, verbose = FALSE){
  res.init.hc <- hc(don,
                    modelName = "VVV",
                    use = "STD")
  res.mclust <- Mclust(don, G = G, modelNames = modelNames, verbose = verbose, initialization = list(hcPairs = res.init.hc))
  res.mclust.boot <- MclustBootstrap(object = res.mclust,nboot=1,type="bs",verbose=verbose)
  res.mclust$param$pro<-res.mclust.boot$pro
  res.mclust$param$mean<-matrix(res.mclust.boot$mean, nrow=dim(res.mclust.boot$mean)[3], byrow=TRUE)
  res.mclust$param$variance$sigma<-array(res.mclust.boot$variance,
                                         dim=dim(res.mclust$param$variance$sigma),
                                         dimnames = dimnames(res.mclust$param$variance$sigma))
  return(res.mclust)
}