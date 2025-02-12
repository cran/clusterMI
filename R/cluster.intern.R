#' Apply clustering method after multiple imputation
#' 
#' @description From a list of imputed datasets \code{clusterMI} performs cluster analysis on each imputed data set.
#' @details Performs cluster analysis (according to the \code{method.clustering} argument). For achieving this goal, the function uses as an input an output from the \code{imputedata} function and applies the cluster analysis method on each imputed data set
#'  
#' Step 1 can be tuned by specifying the cluster analysis method used (\code{method.clustering} argument).
#' If \code{method.clustering = "kmeans"} or \code{"pam"}, then the number of clusters can be specified by tuning the \code{nb.clust} argument. By default, the same number as the one used for imputation is used.
#' The number of random initializations can also be tuned through the \code{nstart.kmeans} argument.
#' If \code{method.clustering = "hclust"} (hierarchical clustering), the method used can be specified (see \code{\link[stats]{hclust}}). By default \code{"average"} is used. Furthermore, the number of clusters can be specified, but it can also be automatically chosen if \code{nb.clust} < 0.
#' If \code{method.clustering = "mixture"} (model-based clustering using gaussian mixture models), the model to be fitted can be tuned by modifying the \code{modelNames} argument (see \code{\link[mclust]{Mclust}}).
#' If \code{method.clustering = "cmeans"} (clustering using the fuzzy c-means algorithm), then the fuzziness parameter can be modfied by tuning the\code{m.cmeans} argument. By default, \code{m.cmeans = 2}.
#' 
#' 
#' Can be performed in parallel by specifying the number of CPU cores (\code{nnodes} argument).
#'
#' @param res.imp a list of imputed data sets
#' @param method.clustering a single string specifying the clustering algorithm used ("kmeans", "pam", "clara", "hclust" or "mixture","cmeans")
#' @param scaling boolean. If TRUE, variables are scaled. Default value is TRUE
#' @param nb.clust an integer specifying the number of clusters
#' @param method.hclust character string defining the clustering method for hierarchical clustering (required only if method.clustering = "hclust")
#' @param method.dist character string defining the method use for computing dissimilarity matrices in hierarchical clustering (required only if method.clustering = "hclust")
#' @param modelNames character string indicating the models to be fitted in the EM phase of clustering (required only if method.clustering = "mixture"). By default modelNames = NULL.
#' @param modelName.hc A character string indicating the model to be used in model-based agglomerative hierarchical clustering.(required only if method.clustering = "mixture"). By default modelNames.hc = "VVV".
#' @param nstart.kmeans how many random sets should be chosen for kmeans initalization. Default value is 100 (required only if method.clustering = "kmeans")
#' @param iter.max.kmeans how many iterations should be chosen for kmeans. Default value is 10 (required only if method.clustering = "kmeans")
#' @param m.cmeans degree of fuzzification in cmeans clustering. By default m.cmeans = 2
#' @param samples.clara number of samples to be drawn from the dataset when performing clustering using clara algorithm. Default value is 500.
#' @param verbose logical
#' @return A list with clustering results
#' @seealso \code{\link[stats]{hclust}}, \code{\link[mclust]{Mclust}}, \code{\link{imputedata}}, \code{\link[e1071]{cmeans}},\code{\link[stats]{dist}}
#' @keywords internal
#' @importFrom e1071 cmeans
#' @importFrom stats kmeans cutree hclust dist
#' @importFrom FactoMineR PCA MCA FAMD
#' @importFrom mclust Mclust hc hcEEE hcVVV hcVII hcEII
#' @importFrom diceR CSPA
#' @importFrom fpc kmeansCBI nselectboot claraCBI noisemclustCBI hclustCBI
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' data(wine)
#' 
#' require(parallel)
#' set.seed(123456)
#' ref <- wine$cult
#' nb.clust <- 3
#' m <- 5 # number of imputed data sets. Should be larger in practice
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- prodna(wine.na)
#' 
#' #imputation
#' res.imp <- imputedata(data.na = wine.na, nb.clust = nb.clust, m = m)
#' \donttest{
#' #analysis by kmeans and pooling
#' nnodes <- 2 # parallel::detectCores()
#' res.pool <- clusterMI(res.imp, nnodes = nnodes)
#' 
#' res.pool$instability
#' table(ref, res.pool$part)
#' }

cluster.intern<-function(res.imp,
                    method.clustering = "kmeans",
                    scaling = TRUE,
                    nb.clust = NULL,
                    method.hclust = "average",
                    method.dist = "euclidean",
                    modelNames = NULL,
                    modelName.hc = "VVV",
                    nstart.kmeans = 100,
                    iter.max.kmeans = 10,
                    m.cmeans = 2,
                    samples.clara=500,
                    verbose = FALSE
                    ){

  if(method.clustering=="kmeans"){
    res.clust<-lapply(res.imp,FUN=function(xx, scaling,nb.clust,nstart.kmeans,iter.max.kmeans){
      if(scaling){
        res.out<-kmeans(scale(xx),centers=nb.clust,nstart = nstart.kmeans,iter.max=iter.max.kmeans)
      }else{
        res.out<-kmeans(xx,centers=nb.clust,nstart = nstart.kmeans,iter.max=iter.max.kmeans)
      }
      return(res.out)
    },scaling=scaling,nb.clust=nb.clust,nstart.kmeans=nstart.kmeans,iter.max.kmeans=iter.max.kmeans)
    res.clust.part<-sapply(res.clust,"[[","cluster")
    res.clust.part.list<-lapply(res.clust,"[[","cluster")
  }
  
  
  if(method.clustering=="cmeans"){
    res.clust<-lapply(res.imp,FUN=function(xx, scaling,nb.clust, m.cmeans){
      if(scaling){
        res.out<-cmeans(scale(xx), centers = nb.clust, m = m.cmeans)
      }else{
        res.out<-cmeans(xx,centers=nb.clust, m = m.cmeans)
      }
      return(res.out)
    },scaling=scaling,nb.clust=nb.clust, m.cmeans=m.cmeans)
    res.clust.part<-sapply(res.clust,"[[","cluster")
    res.clust.part.list<-lapply(res.clust,"[[","cluster")
  }
  
  if(method.clustering=="clara"){
    if(scaling){
      res.clust<-lapply(res.imp,FUN=claraCBI,k=nb.clust,usepam=FALSE,diss=FALSE,stand=TRUE,samples=samples.clara)
    }else if (!scaling){
      res.clust<-lapply(res.imp,FUN=claraCBI,k=nb.clust,usepam=FALSE,diss=FALSE,stand=FALSE,samples=samples.clara)
    }
    res.clust.part.list<-lapply(res.clust,FUN = function(xx){xx[["result"]][["clustering"]]})
    res.clust.part<-do.call(cbind,res.clust.part.list)
  }
  if(method.clustering=="pam"){
    if(scaling){
      res.clust<-lapply(res.imp,FUN=claraCBI,k=nb.clust,usepam=TRUE,diss=FALSE,stand=TRUE)
    }else if (!scaling){
      res.clust<-lapply(res.imp,FUN=claraCBI,k=nb.clust,usepam=TRUE,diss=FALSE,stand=FALSE)
    }
    res.clust.part.list<-lapply(res.clust,FUN = function(xx){xx[["result"]][["clustering"]]})
    res.clust.part<-do.call(cbind,res.clust.part.list)
  }
  
  if(method.clustering=="hclust"){
    res.clust.part.list<-lapply(res.imp,
                                FUN=function(xx,scaling,method.dist,method.hclust,nb.clust){
                                  if(scaling){xxs<-scale(xx,center=FALSE)
                                  } else if (!scaling){
                                    xxs<-xx
                                  }
                                  d <-dist(xxs,method=method.dist)
                                  
                                  res.hclust<- hclust(d=d,method=method.hclust)
                                  
                                  if(nb.clust>0){
                                    part <- cutree(res.hclust,k=nb.clust)
                                  }else if (nb.clust <0){
                                    kgrid<-seq.int(from = 2,to = min(nrow(xx)/10,15))
                                    names(kgrid)<-kgrid
                                    Sil <- sapply(kgrid,FUN = function(k,d){
                                      Silhouette.intern(d = d,cl = cutree(res.hclust,k=k))
                                    },d=d)
                                    names(Sil)<-kgrid
                                    kopt <- as.numeric(names(which.max(Sil)))
                                    part <- cutree(res.hclust,k=kopt)
                                  }
                                  return(part)
                                },method.dist=method.dist,method.hclust= method.hclust,scaling=scaling,nb.clust=nb.clust)
    res.clust.part<-do.call(cbind,res.clust.part.list)
  }
  
  if(method.clustering=="mixture"){
    if(verbose){cat("\n");pb <- utils::txtProgressBar(style = 3)}else{pb<-NULL}
    res.clust <- mapply(res.imp, seq(length(res.imp)), 
                        FUN = function(xx, indice, scaling, nb.clust, 
                                       modelNames, modelName, pb, indicemax, seed=NULL) {
                          
                          if (scaling) {
                            yy <- scale(xx)
                          }
                          else {
                            yy <- xx
                          }
                          if(is.null(seed)){
                            partition<-suppressWarnings({
                              kmeans(yy, centers = min(nrow(yy), 100), 
                                     nstart = 100)$cluster
                            })
                          }else{
                            partition<-suppressWarnings({
                              with_seed(seed=seed,
                                        kmeans(yy, centers = min(nrow(yy), 100), 
                                               nstart = 100)$cluster)
                            })
                          }
                          
                          res.init.hc <- hc(yy, modelName = modelName, 
                                            use = "VARS", partition = partition)
                          res.out <- Mclust(data = yy, G = nb.clust, 
                                            modelNames = modelNames, verbose = FALSE, 
                                            initialization = list(hcPairs = res.init.hc))
                          try(setTxtProgressBar(pb, indice/indicemax),silent=TRUE)
                          return(res.out)
                        }, MoreArgs = list(scaling = scaling, 
                                           nb.clust = nb.clust, modelNames = modelNames, 
                                           modelName = modelName.hc, pb = pb, indicemax = length(res.imp), seed=NULL), 
                        SIMPLIFY = FALSE)
    if(verbose){close.txtProgressBar(pb)}
    
    res.clust.part<-sapply(res.clust,"[[","classification")
    res.clust.part.list<-lapply(res.clust,"[[","classification")
  }
    res.out<-list(part.list=res.clust.part.list,part.matrix=res.clust.part)

  return(res.out)
}