#' initialize fastnmf
#'
#' @param method a vector giving initialisation methods among "BOK", "kmeans", "minibatchkmeans" "sample"
#' @param Mtilde a fuzzy connectivity matrix
#' @param listpart a list of partitions (required for kmeans and minibatchkmeans method)
#' @param parameter.kmeans list of arguments for kmeans function
#' @param parameter.minibatchkmeans list of arguments for MiniBatchKmeans function
#' @importFrom ClusterR MiniBatchKmeans predict_MBatchKMeans
#' @importFrom Rfast Tcrossprod
#' @importFrom stats kmeans
#' @keywords internal
#' 
initfastnmf <- function(method,Mtilde,listpart=NULL,nb.clust=NULL,
                        parameter.kmeans=list(nstart=100,
                                              iter.max=50,
                                              algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                                                            "MacQueen"),
                                              trace = FALSE),
                        parameter.minibatchkmeans=list(batch_size = 10,
                                                       num_init = 1,
                                                       max_iters = 50, 
                                                       init_fraction = 1,
                                                       initializer = "kmeans++",
                                                       early_stop_iter = 10, 
                                                       verbose = FALSE,
                                                       CENTROIDS = NULL,
                                                       tol = 1e-04,
                                                       tol_optimal_init = 0.3, 
                                                       seed = 1)){
  if(method=="BOK"){
    if(is.null(listpart)){stop("method is BOK, but listpart is NULL")}else if (!is.null(listpart)){
      res.out <- listpart[[which.min(sapply(listpart,FUN=function(part,Mtilde){
        H<-class2ind(factor(part))
        res.out <- sum((Mtilde-Tcrossprod(H, H))^2)
        return(res.out)
      },Mtilde=Mtilde)
      )]]
    }
  }else if(method=="kmeans"){
    res.kmeans <- try(kmeans(Mtilde,
                             centers = nb.clust, 
                             nstart = parameter.kmeans$nstart,
                             iter.max = parameter.kmeans$iter.max,
                             algorithm = parameter.kmeans$algorithm,
                             trace=parameter.kmeans$trace))
    if (!inherits(res.kmeans,"try-error")){
      res.out<- res.kmeans$cluster
    }else{
      res.out <- NULL
    }
  }else if(method=="MiniBatchKmeans"){
    res.kmeans<-try(MiniBatchKmeans(Mtilde,
                                    clusters = nb.clust,
                                    batch_size = parameter.minibatchkmeans$batch_size,
                                    num_init = parameter.minibatchkmeans$nstart,
                                    max_iters = parameter.minibatchkmeans$iter.max,
                                    init_fraction = parameter.minibatchkmeans$init_fraction,
                                    initializer = parameter.minibatchkmeans$initializer,
                                    early_stop_iter = parameter.minibatchkmeans$early_stop_iter,
                                    verbose = parameter.minibatchkmeans$verbose,
                                    CENTROIDS = parameter.minibatchkmeans$CENTROIDS,
                                    tol=parameter.minibatchkmeans$tol,
                                    tol_optimal_init = parameter.minibatchkmeans$tol_optimal_init,
                                    seed = parameter.minibatchkmeans$seed))
    if(!inherits(res.kmeans,"try-error")){
      res.out<- predict_MBatchKMeans(Mtilde,res.kmeans$centroids, fuzzy=FALSE)
    }
  }else if(method=="sample"){
    res.out <- sample(seq(nb.clust),size = ncol(Mtilde),replace=TRUE)
  }else{stop("initialization method is unknown")}
}
