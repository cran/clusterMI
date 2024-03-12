#' Consensus clustering using non-negative matrix factorization
#' 
#' @description From a list of partitions \code{fastnmf} pools partition as proposed in Li and Ding (2007) <doi:10.1109/ICDM.2007.98>.
#' @return a list of 5 objets
#'  \item{Htilde}{A fuzzy disjunctive table}
#'  \item{S}{A positive matrix}
#'  \item{Mtilde}{The average of connectivity matrices}
#'  \item{crit}{A vector with the optimized criterion at each iteration}
#'  \item{cluster}{the consensus partition in nb.clust classes}
#' @param listpart a list of partitions
#' @param nb.clust an integer specifying the number of clusters
#' @param threshold a real specifying when the NMF algorithm is stoped. Default value is 10^(-5)
#' @param printflag a boolean. If TRUE, nmf will print messages on console. Default value is TRUE
#' @param nstart how many random sets should be chosen for kmeans initalization. Default value is 100
#' @param early_stop_iter continue that many iterations after calculation of the best within-cluster-sum-of-squared-error. Default value is 10. See MiniBatchKmeans help page.
#' @param initializer the method of initialization. One of, optimal_init, quantile_init, kmeans++ and random. See MiniBatchKmeans help page.
#' @param batch_size the size of the mini batches for kmeans clustering. Default value is NULL.
#' @param iter.max the maximum number of iterations allowed for kmeans. Default value is 50
#' @seealso \code{\link[stats]{kmeans}} \code{\link[ClusterR]{MiniBatchKmeans}}
#' @details fastnmf performs consensus clustering using non-negative matrix factorization following Li and Ding (2007) <doi:10.1109/ICDM.2007.98>. The set of partitions that are aggregated needs to be given as a list where each element is a vector of numeric values. Note that the number of classes for each partition can vary. The number of classes for the consensus partition should be given using the nb.clust argument. The NMF algorithm is iterative and required an initial partition. This latter is based on kmeans clustering on the average of connectivity matrices. If batchsize is NULL, then kmeans clustering is performed using \code{nstart} initial values and \code{iter.max} iterations. Otherwise, Mini Batch Kmeans is used. This algorithm could be faster than kmeans if the number of invididuals is large.
#' @references T. Li, C. Ding, and M. I. Jordan (2007) Solving consensus and semi-supervised clustering problems using nonnegative matrix factorization.  In Proceedings of the 2007 Seventh IEEE International Conference on Data Mining, ICDM'07, page 577-582, USA. IEEE Computer Society. <doi:10.1109/ICDM.2007.98>
#' @importFrom ClusterR MiniBatchKmeans predict_MBatchKMeans
#' @importFrom Rfast Crossprod Tcrossprod mat.mult
#' @export
#' @examples
#' data(wine)
#' require(clustrd)
#' set.seed(123456)
#' ref <- wine$cult
#' nb.clust <- 3
#' m <- 3 # number of imputed data sets. Should be larger in practice
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- as.matrix(wine.na)
#' wine.na[sample(seq(length(wine.na)), size = ceiling(length(wine.na)/3))] <- NA
#' 
#' #imputation
#' res.imp <- imputedata(data.na = wine.na, nb.clust = nb.clust, m = m)
#' 
#' #analysis using reduced kmeans
#'
#' ## apply the cluspca function on each imputed data set
#' res.ana.rkm <- lapply(res.imp$res.imp,
#'                       FUN = cluspca,
#'                       nclus = nb.clust,
#'                       ndim = 2,
#'                       method= "RKM")
#' ## extract the set of partitions (under "list" format)
#' res.ana.rkm <-lapply(res.ana.rkm,"[[","cluster")
#' 
#' # pooling by NMF
#' res.pool.rkm <- fastnmf(res.ana.rkm, nb.clust = nb.clust)$clust
#' 


fastnmf<-function(listpart,nb.clust,threshold=10^(-5),printflag=TRUE,nstart=100,early_stop_iter=10,initializer = 'random',batch_size = NULL,iter.max=50){
  fcritsave<-function(Mtilde,Htilde,S){
    tcSHtilde<-Tcrossprod(S,Htilde)
    tmp<- Mtilde-mat.mult(Htilde,tcSHtilde)
    res<-sqrt(sum(apply(tmp,2,crossprod)))
    return(res)
  }
  
  listpart.intern<-listpart[!sapply(listpart,is.null)]
  if(length(listpart.intern)==1){
    res<-list(Htilde=NULL,S=NULL,crit=NULL,Mtilde=NULL,cluster=listpart.intern[[1]])
    return(res)
  }else if(length(listpart.intern)==0){ res<-list(Htilde=NULL,S=NULL,crit=NULL,Mtilde=NULL,cluster=NULL);return(res)}
  if(is.null(batch_size)&(length(listpart.intern[[1]])>=1000)){warning("The numbers of individuals is large. You should specify the batch_size argument to accelerate the initialisation step.")}
  Mtilde<-matrix(0,length(listpart.intern[[1]]),length(listpart.intern[[1]]))#matrice d'adjacence rempli de 0
  for(ii in seq(length(listpart.intern))){
    Mtilde<-Mtilde+tcrossprod(class2ind(factor(listpart.intern[[ii]])))#somme de toutes les matrices d'adjacences
  }
  Mtilde<-Mtilde/length(listpart.intern)#moyenne des matrices d'adjacences
  
  #initialisation
  if(!is.null(batch_size)){
    #minibatch
    res.kmeans<-try(MiniBatchKmeans(Mtilde, clusters = nb.clust, batch_size = batch_size, num_init = nstart, max_iters = iter.max, 
                                    initializer = initializer, early_stop_iter = early_stop_iter,
                                    verbose = FALSE))
    if("try-error"%in% class(res.kmeans)){
      res.kmeans<-sample(seq(nb.clust),size = ncol(Mtilde),replace=TRUE)
    }else{
      res.kmeans<- predict_MBatchKMeans(Mtilde,res.kmeans$centroids, fuzzy=FALSE)
    }
  }else{
    #kmeans (stats R package)
    res.kmeans<-try(kmeans(Mtilde, centers = nb.clust,  nstart = nstart, iter.max = iter.max))
    if("try-error"%in% class(res.kmeans)){
      res.kmeans<-sample(seq(nb.clust),size = ncol(Mtilde),replace=TRUE)
    }else{
      res.kmeans<- res.kmeans$cluster
    }
  }
  H<-class2ind(factor(res.kmeans))
  S<-Crossprod(H,H)
  Htilde<-mat.mult(H,diag(diag(1/sqrt(S)),nb.clust,nb.clust))
  continue<-TRUE
  critsave<-fcritsave(Mtilde = Mtilde,Htilde = Htilde,S = S)
  comp<-1
  while(continue){
    if(printflag){cat(comp,"...")}
    #Htilde update
    MHS <- mat.mult(Mtilde,Htilde)%*%S
    cpHtilde <-Tcrossprod(Htilde,Htilde)
    multHtilde<-sqrt(MHS/(mat.mult(cpHtilde,MHS)))
    multHtilde[is.nan(multHtilde)]<-0
    Htilde<-Htilde*multHtilde
    #S update
    cpHtilde <- Crossprod(Htilde,Htilde)
    cpHtildeMtilde<-Crossprod(Htilde,Mtilde)
    multS<-sqrt((mat.mult(cpHtildeMtilde,Htilde))/(cpHtilde%*%S%*%cpHtilde))
    multS[is.nan(multS)]<-0
    multS[is.infinite(multS)]<-0
    S<-S*multS
    
    critsave<-c(critsave,fcritsave(Mtilde = Mtilde,Htilde = Htilde,S = S))
    comp<-comp+1
    diffcrit<-(critsave[comp-1]-critsave[comp])/critsave[comp-1]
    continue<-(diffcrit>=threshold)
  }
  if(printflag){cat("done! \n")}
  
  res<-list(Htilde=Htilde,S=S,Mtilde=Mtilde,crit=critsave,cluster=apply(Htilde,1,which.max))
  return(res)
}