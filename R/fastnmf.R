#' Consensus clustering using non-negative matrix factorization
#' 
#' @description From a list of partitions \code{fastnmf} pools partition as proposed in Li and Ding (2007) <doi:10.1109/ICDM.2007.98>.
#' @return For each initialisation method, a list of 5 objets is returned
#'  \item{Htilde}{A fuzzy disjunctive table}
#'  \item{S}{A positive matrix}
#'  \item{Mtilde}{The average of connectivity matrices}
#'  \item{crit}{A vector with the optimized criterion at each iteration}
#'  \item{cluster}{the consensus partition in nb.clust classes}
#'  In addition, the best initialisation method is returned
#' @param listpart a list of partitions
#' @param nb.clust an integer specifying the number of clusters
#' @param method.init a vector giving initialisation methods used among "BOK", "kmeans", "minibatchkmeans" "sample". See details.
#' @param threshold a real specifying when the NMF algorithm is stoped. Default value is 10^(-5)
#' @param printflag a boolean. If TRUE, nmf will print messages on console. Default value is TRUE
#' @param parameter.kmeans a list of arguments for kmeans function. See keans help page. 
#' @param parameter.minibatchkmeans list of arguments for MiniBatchKmeans function. See MiniBatchKmeans help page.
#' @seealso \code{\link[stats]{kmeans}} \code{\link[ClusterR]{MiniBatchKmeans}}
#' @details fastnmf performs consensus clustering using non-negative matrix factorization following Li and Ding (2007) <doi:10.1109/ICDM.2007.98>. The set of partitions that are aggregated needs to be given as a list where each element is a vector of numeric values. Note that the number of classes for each partition can vary. The number of classes for the consensus partition should be given using the \code{nb.clust} argument. The NMF algorithm is iterative and required an initial partition. This latter is specified by \code{method.init}. \code{method.init="BOK"} means the partition considered is a partition from \code{listpart} which minimizes the NMF criterion. Alternative methods are "kmeans", "minibathckmeans" or "sample". If \code{method.init} = "kmeans" (or "minibatchkmeans"), then clustering on the average of connectivity matrices is performed by kmeans (or "minibatchkmeans"). Mini Batch Kmeans could be faster than kmeans if the number of invididuals is large. If \code{method.init} = "sample", then a random partition is drawn. If \code{method.init} is a vector of several characters, then several initialization methods are considered and the best method is returned. By default, \code{method.init= c("BOK", "kmeans")}.
#' @references T. Li, C. Ding, and M. I. Jordan (2007) Solving consensus and semi-supervised clustering problems using nonnegative matrix factorization.  In Proceedings of the 2007 Seventh IEEE International Conference on Data Mining, ICDM'07, page 577-582, USA. IEEE Computer Society. <doi:10.1109/ICDM.2007.98>
#' @importFrom ClusterR MiniBatchKmeans predict_MBatchKMeans
#' @importFrom Rfast Crossprod Tcrossprod mat.mult
#' @importFrom stats kmeans
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
#' wine.na <- prodna(wine.na)
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
#' res.pool.rkm <- fastnmf(res.ana.rkm, nb.clust = nb.clust)
#' ## extract the partition corresponding to the best initialisation
#' part <- res.pool.rkm$best$clust
#' 


fastnmf<-function(listpart,
                  nb.clust,
                  method.init=c("BOK","kmeans"),
                  threshold=10^(-5),
                  printflag=TRUE,
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
  fcritsave<-function(Mtilde,Htilde,S){
    tcSHtilde<-Tcrossprod(S,Htilde)
    tmp<- Mtilde-mat.mult(Htilde,tcSHtilde)
    res<-sqrt(sum(apply(tmp,2,crossprod)))
    return(res)
  }
  
  listpart.intern<-listpart[!sapply(listpart,is.null)]
  if(length(listpart.intern)==1){
    res<-list(Htilde=NULL,S=NULL,crit=NULL,Mtilde=NULL,cluster=listpart.intern[[1]])
    res$best <- res
    return(res)
  }else if(length(listpart.intern)==0){ res<-list(Htilde=NULL,S=NULL,crit=NULL,Mtilde=NULL,cluster=NULL);return(res)}
  if((!("MiniBatchKmeans"%in%method.init))&("kmeans"%in%method.init)&(length(listpart.intern[[1]])>=1000)){warning("The numbers of individuals is large. You should specify the batch_size argument to accelerate the initialisation step.")}
  Mtilde<-matrix(0,length(listpart.intern[[1]]),length(listpart.intern[[1]]))#matrice d'adjacence rempli de 0
  for(ii in seq(length(listpart.intern))){
    Mtilde<-Mtilde+tcrossprod(class2ind(factor(listpart.intern[[ii]])))#somme de toutes les matrices d'adjacences
  }
  Mtilde<-Mtilde/length(listpart.intern)#moyenne des matrices d'adjacences
  
  #initialisation
  res.kmeans <- sapply(method.init,
                       initfastnmf,
                       Mtilde = Mtilde,
                       nb.clust = nb.clust,
                       listpart = listpart,
                       parameter.kmeans=parameter.kmeans,
                       parameter.minibatchkmeans=parameter.minibatchkmeans,
                       simplify=FALSE,USE.NAMES = TRUE)
  
  res.nmf <- sapply(method.init,FUN=function(method.init,initparts,nb.clust,threshold,printflag){
    if(printflag){cat("Initialisation by", method.init, " ")}
    initpart <- initparts[[method.init]]
    H<-class2ind(factor(initpart))
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
  },initparts = res.kmeans,nb.clust=nb.clust,threshold=threshold,printflag=printflag,simplify=FALSE)
  
  critbyinit <- sapply(res.nmf,"[[","crit", simplify = FALSE)
  
  res.nmf$best <- res.nmf[[which.min(sapply(critbyinit,min,simplify=TRUE))]]
  return(res.nmf)
}