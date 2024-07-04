#' Cluster analysis and pooling after multiple imputation
#' 
#' @description From a list of imputed datasets \code{clusterMI} performs cluster analysis on each imputed data set, estimates the instability of each partition using bootstrap (following Fang, Y. and Wang, J., 2012 <doi:10.1016/j.csda.2011.09.003>) and pools results as proposed in Audigier and Niang (2022) <doi:10.1007/s11634-022-00519-1>.
#' @details \code{clusterMI} performs cluster analysis (according to the \code{method.clustering} argument) and pooling after multiple imputation. For achieving this goal, the \code{clusterMI} function uses as an input an output from the \code{imputedata} function and then
#' \enumerate{
#'  \item applies the cluster analysis method on each imputed data set
#'  \item pools contributory partitions using non-negative matrix factorization
#'  \item computes the instability of each partition by bootstrap
#'  \item computes the total instability
#'  }
#'  
#' Step 1 can be tuned by specifying the cluster analysis method used (\code{method.clustering} argument).
#' If \code{method.clustering = "kmeans"} or \code{"pam"}, then the number of clusters can be specified by tuning the \code{nb.clust} argument. By default, the same number as the one used for imputation is used.
#' The number of random initializations can also be tuned through the \code{nstart.kmeans} argument.
#' If \code{method.clustering = "hclust"} (hierarchical clustering), the method used can be specified (see \code{\link[stats]{hclust}}). By default \code{"average"} is used. Furthermore, the number of clusters can be specified, but it can also be automatically chosen if \code{nb.clust} < 0.
#' If \code{method.clustering = "mixture"} (model-based clustering using gaussian mixture models), the model to be fitted can be tuned by modifying the \code{modelNames} argument (see \code{\link[mclust]{Mclust}}).
#' If \code{method.clustering = "cmeans"} (clustering using the fuzzy c-means algorithm), then the fuzziness parameter can be modfied by tuning the\code{m.cmeans} argument. By default, \code{m.cmeans = 2}.
#' 
#' Step 2 performs consensus clustering by Non-Negative Matrix Factorization, following Li and Ding (2007) <doi:10.1109/ICDM.2007.98>.
#' 
#' Step 3 applies the \code{\link[fpc]{nselectboot}} function on each imputed data set and returns the instability of each cluster obtained at step 1. The method is based on bootstrap sampling, followong Fang, Y. and Wang, J. (2012) <doi:10.1016/j.csda.2011.09.003>. The number of iterations can be tuned using the \code{Cboot} argument.
#' 
#' Step 4 averages the previous instability measures given a within instability (\code{Ubar}), computes a between instability (\code{B}) and a total instability (\code{T} = B + Ubar). See Audigier and Niang (2022) <doi:10.1007/s11634-022-00519-1> for details.
#' 
#' All steps can be performed in parallel by specifying the number of CPU cores (\code{nnodes} argument). Steps 3 and 4 are more time consuming. To compute only steps 1 and 2 use \code{instability = FALSE}.
#'
#' @param output an output from the imputedata function
#' @param method.clustering a single string specifying the clustering algorithm used ("kmeans", "pam", "clara", "hclust" or "mixture","cmeans")
#' @param method.consensus a single string specifying the consensus method used to pool the contributory partitions ("NMF" or "CSPA")
#' @param scaling boolean. If TRUE, variables are scaled. Default value is TRUE
#' @param nb.clust an integer specifying the number of clusters
#' @param Cboot an integer specifying the number of bootstrap replications. Default value is 50
#' @param method.hclust character string defining the clustering method for hierarchical clustering (required only if method.clustering = "hclust")
#' @param method.dist character string defining the method use for computing dissimilarity matrices in hierarchical clustering (required only if method.clustering = "hclust")
#' @param modelNames character string indicating the models to be fitted in the EM phase of clustering (required only if method.clustering = "mixture"). By default modelNames = NULL.
#' @param modelName.hc A character string indicating the model to be used in model-based agglomerative hierarchical clustering.(required only if method.clustering = "mixture"). By default modelNames.hc = "VVV".
#' @param nstart.kmeans how many random sets should be chosen for kmeans initalization. Default value is 100 (required only if method.clustering = "kmeans")
#' @param iter.max.kmeans how many iterations should be chosen for kmeans. Default value is 10 (required only if method.clustering = "kmeans")
#' @param m.cmeans degree of fuzzification in cmeans clustering. By default m.cmeans = 2
#' @param samples.clara number of samples to be drawn from the dataset when performing clustering using clara algorithm. Default value is 500.
#' @param nnodes number of CPU cores for parallel computing. By default, nnodes = 1
#' @param instability a boolean indicating if cluster instability must be computed. Default value is TRUE
#' @param verbose a boolean. If TRUE, a message is printed at each step. Default value is TRUE
#' @param nmf.threshold Default value is 10^(-5),
#' @param nmf.nstart Default value is 100,
#' @param nmf.early_stop_iter Default value is 10,
#' @param nmf.initializer Default value is 'random',
#' @param nmf.batch_size Default value is 20,
#' @param nmf.iter.max Default value is 50
#' @return A list with three objects
#' \item{part}{the consensus partition}
#' \item{instability}{a list of four objects: \code{U} the within instability measure for each imputed data set, \code{Ubar} the associated average, \code{B} the between instability measure, \code{Tot} the total instability measure}
#' \item{call}{the matching call}
#' @seealso \code{\link[stats]{hclust}}, \code{\link[fpc]{nselectboot}}, \code{\link[mclust]{Mclust}}, \code{\link{imputedata}}, \code{\link[e1071]{cmeans}}
#' @aliases clusterMI
#' @references 
#'   Audigier, V. and Niang, N. (2022) Clustering with missing data: which equivalent for Rubin's rules? Advances in Data Analysis and Classification <doi:10.1007/s11634-022-00519-1>
#'   
#'   Fang, Y. and Wang, J. (2012) Selection of the number of clusters via the bootstrap method. Computational Statistics and Data Analysis, 56, 468-477. <doi:10.1016/j.csda.2011.09.003>
#'   
#'   T. Li, C. Ding, and M. I. Jordan (2007) Solving consensus and semi-supervised clustering problems using nonnegative matrix factorization.  In Proceedings of the 2007 Seventh IEEE International Conference on Data Mining, ICDM'07, page 577-582, USA. IEEE Computer Society. <doi:10.1109/ICDM.2007.98>
#' @export
#' @importFrom FactoMineR PCA MCA FAMD
#' @importFrom diceR CSPA
#' @importFrom fpc kmeansCBI nselectboot claraCBI noisemclustCBI hclustCBI
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

clusterMI<-function(output,
                  method.clustering = "kmeans",
                  method.consensus = "NMF",
                  scaling = TRUE,
                  nb.clust = NULL,
                  Cboot = 50,
                  method.hclust = "average",
                  method.dist = "euclidean",
                  modelNames = NULL,
                  modelName.hc = "VVV",
                  nstart.kmeans = 100,
                  iter.max.kmeans = 10,
                  m.cmeans = 2,
                  samples.clara=500,
                  nnodes = 1,
                  instability = TRUE,
                  verbose = TRUE,
                  nmf.threshold=10^(-5),
                  nmf.nstart=100,
                  nmf.early_stop_iter=10,
                  nmf.initializer = 'random',
                  nmf.batch_size = NULL,
                  nmf.iter.max=50){
  
  res.imp<-output$res.imp
  if(any(sapply(lapply(res.imp,is.na),sum)>0)){
    warning("NA in imputed datasets. Incomplete variables have been removed.")
    res.imp<-lapply(res.imp,FUN=function(xx){
      is.complete<-apply(!is.na(xx),2,all)
      xx.out<-xx[,is.complete,drop=FALSE]
      return(xx.out)})
    }
  if(is.null(nb.clust)){nb.clust.intern<-output$call$nb.clust}else{nb.clust.intern<-nb.clust}
  if(any(sapply(res.imp,FUN=nrow)<=15)|("try-error"%in%class(res.imp))){
    res.out<-list(part=NA,
                  instability=list(Ubar=NA,U=NA,B=NA),
                  call=list(output=output,
                            method.consensus=method.consensus,
                            method.clustering=method.clustering,
                            scaling=scaling,
                            nb.clust=nb.clust,
                            Cboot=Cboot,
                            method.hclust=method.hclust,
                            method.dist=method.dist,
                            modelNames=modelNames,
                            modelName.hc=modelName.hc,
                            nstart.kmeans=nstart.kmeans,
                            iter.max.kmeans=iter.max.kmeans,
                            m.cmeans = m.cmeans,
                            nnodes=nnodes,
                            instability=instability,
                            nmf.threshold=nmf.threshold,
                            nmf.nstart=nmf.nstart,
                            nmf.early_stop_iter = nmf.early_stop_iter,
                            nmf.initializer =  nmf.initializer,
                            nmf.batch_size = nmf.batch_size,
                            nmf.iter.max = nmf.iter.max
                            ))
    return(res.out)
  }
  if(output$call$data.type=="categorical"){
    warning("Data are categorical, cluster analysis is performed on principal components")
    if(scaling){warning("scaling = TRUE does not make sense with categorical data")}
    res.imp.intern<-lapply(res.imp,FUN=function(xx){as.data.frame(MCA(xx,ncp=Inf,graph=FALSE)$ind$coord)})
    scaling.intern<-FALSE
  }else if(output$call$data.type=="mixed"){
    warning("Data are mixed, cluster analysis is performed on principal components")
    if(scaling){warning("scaling = TRUE does not make sense with mixed data")}
    
    res.imp.intern<-lapply(res.imp,FUN=function(xx){
      if(all(sapply(xx,is.factor)|sapply(xx,is.ordered))){
        # il ne reste que des quali après suppression des variables incomplètes
        res.out <- as.data.frame(MCA(xx,ncp=Inf,graph=FALSE)$ind$coord)
      } else if(all(sapply(xx,is.numeric)|sapply(xx,is.integer))){
        # il ne reste que des quanti après suppression des variables incomplètes
        res.out <- as.data.frame(PCA(xx,ncp=Inf,graph=FALSE)$ind$coord)
          }else{
      res.out <- as.data.frame(FAMD(xx,ncp=Inf,graph=FALSE)$ind$coord)
          }
      return(res.out)
      })
    scaling.intern<-FALSE
  }else if(output$call$data.type=="continuous"){
    res.imp.intern<-res.imp
    scaling.intern<-scaling
        }

    #Analysis
   if(verbose){cat("Cluster analysis (by ", method.clustering,")...",sep="")}
  
  res.ana <- cluster.intern(res.imp = res.imp.intern,
                            method.clustering = method.clustering,
                            scaling = scaling.intern,
                            nb.clust = nb.clust.intern,
                            method.hclust = method.hclust,
                            method.dist = method.dist,
                            modelNames = modelNames,
                            modelName.hc = modelName.hc,
                            nstart.kmeans = nstart.kmeans,
                            iter.max.kmeans = iter.max.kmeans,
                            m.cmeans = m.cmeans,
                            samples.clara=samples.clara,
                            verbose = verbose)
  
  res.clust.part <- res.ana$part.matrix
  res.clust.part.list <- res.ana$part.list
  rm("res.ana")
  
  if(verbose){cat("done!\n")}
  
  if(verbose){cat("Consensus clustering (by ", method.consensus,")...",sep="")}
  

  #Pooling
  
  #partition
  
  if(method.consensus=="CSPA"){
    Dimnames<-list(seq(nrow(res.clust.part)), "R1", paste0("m=",seq(length(res.clust.part.list))), nb.clust.intern)
    
    res.pool<-array(NA,dim=sapply(Dimnames,length),dimnames = Dimnames)
    res.pool[,1,,1]<-res.clust.part
    res.pool<-CSPA(res.pool,k=nb.clust.intern)
  }else if(method.consensus=="NMF"){
    # print(k);print(summary(res.clust.part))
    res.pool<-fastnmf(res.clust.part.list,
                      nb.clust=nb.clust.intern,
                      printflag = FALSE,
                      threshold=nmf.threshold,
                      nstart=nmf.nstart,
                      early_stop_iter=nmf.early_stop_iter,
                      initializer = nmf.initializer,
                      batch_size = nmf.batch_size,
                      iter.max=nmf.iter.max)$clust
    
    
  }else{stop("consensus method unknown")}
  
  
  if(verbose){cat("done!\n")}
  
  if(instability){
    if(verbose){cat("Compute within instability...")}
    
    #variance
    # if(length(res.imp)>1){
    if(nnodes>1){
      cl <- parallel::makeCluster(nnodes, type = "PSOCK")
      # parallel::clusterEvalQ(cl, library("fpc"))
      # parallel::clusterEvalQ(cl, library("FactoMineR"))
      
      # clusterEvalQ(cl, library("micemd"))
      
      parallel::clusterExport(cl, list("res.imp.intern","calculintra.intern","method.clustering","Cboot","nb.clust.intern","scaling.intern","method.hclust","method.dist","MCA", "FAMD",
                                       "cmeansCBI.intern","kmeansCBI", "nselectboot", "claraCBI", "noisemclustCBI", "hclustCBI","nstart.kmeans","modelNames","m.cmeans"), envir = environment())
      U<-parallel::parSapply(cl,res.imp.intern,
                             FUN=calculintra.intern,
                             method.clustering=method.clustering,
                             Cboot=Cboot,
                             nb.clust=nb.clust.intern,
                             scaling=scaling.intern,
                             method.agnes=method.hclust,
                             method.dist=method.dist,
                             nstart.kmeans=nstart.kmeans,
                             modelNames=modelNames,
                             m.cmeans=m.cmeans)
      parallel::stopCluster(cl)
      
    }else{
      U<-sapply(res.imp.intern,
                FUN=calculintra.intern,
                method.clustering=method.clustering,
                Cboot=Cboot,nb.clust=nb.clust.intern,
                scaling=scaling.intern,
                method.agnes=method.hclust,
                method.dist = method.dist,
                nstart.kmeans=nstart.kmeans,
                modelNames=modelNames,
                m.cmeans=m.cmeans)
      }
    
    Ubar<-mean(U)
    if(verbose){cat("done!\n")}
    if(verbose){cat("Compute between instability...")}
    if(!all(apply(res.clust.part,2,is.integer))){
      res.clust.part<-apply(res.clust.part,2,FUN=function(xx){as.integer(as.numeric(as.character(xx)))})
    }
    if(length(res.imp)>1){
      if(ncol(res.clust.part)==length(res.imp)){
        B <- calcul_inter(res.clust.part)
      }else{  
        B<-calcul_inter(t(res.clust.part))}
    }else{B<-0}
    
    if(verbose){cat("done!\n")}
    
    res.out<-list(part=res.pool,
                  instability=list(U=U,Ubar=Ubar,B=B, Tot=Ubar+B),
                  call=list(output=output,
                            method.consensus=method.consensus,
                            method.clustering=method.clustering,
                            scaling=scaling,
                            nb.clust=nb.clust.intern,
                            Cboot=Cboot,
                            method.hclust=method.hclust,
                            method.dist= method.dist,
                            modelNames=modelNames,
                            modelName.hc=modelName.hc,
                            nstart.kmeans=nstart.kmeans,
                            iter.max.kmeans=iter.max.kmeans,
                            nnodes=nnodes,
                            instability=instability,
                            m.cmeans=m.cmeans,
                            nmf.threshold= nmf.threshold,
                            nmf.nstart=nmf.nstart,
                            nmf.early_stop_iter=nmf.early_stop_iter,
                            nmf.initializer = nmf.initializer,
                            nmf.batch_size = nmf.batch_size,
                            nmf.iter.max=nmf.iter.max,
                            res.analyse = res.clust.part.list))
  }else{
    res.out<-list(part=res.pool,
                  instability=list(Ubar=NA,U=NA,B=NA),
                  call=list(output=output,
                            method.consensus=method.consensus,
                            method.clustering=method.clustering,
                            scaling=scaling,
                            nb.clust=nb.clust,
                            Cboot=Cboot,
                            method.hclust=method.hclust,
                            method.dist= method.dist,
                            modelNames=modelNames,
                            modelName.hc=modelName.hc,
                            nstart.kmeans=nstart.kmeans,
                            iter.max.kmeans=iter.max.kmeans,
                            nnodes=nnodes,
                            instability=instability,
                            m.cmeans=m.cmeans,
                            nmf.threshold= nmf.threshold,
                            nmf.nstart=nmf.nstart,
                            nmf.early_stop_iter=nmf.early_stop_iter,
                            nmf.initializer = nmf.initializer,
                            nmf.batch_size = nmf.batch_size,
                            nmf.iter.max=nmf.iter.max,
                            res.analyse = res.clust.part.list))
  }
  return(res.out)
}