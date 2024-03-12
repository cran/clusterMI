#' Tune the number of clusters according to the partition instability
#' 
#' @description 
#' \code{choosenbclust} reports the cluster instability according to the number of clusters chosen.
#' 
#' @details 
#'  The \code{choosenbclust} function browses a grid of values for the number of clusters and for each one imputes the data and computes the instability.
#'  
#' @return a list of two objects
#'  \item{nb.clust}{the number of clusters in \code{grid} minimizing the instability}
#'  \item{crit}{a vector indicating the instability for each value in the grid}
#'  
#' @param output an output from the clusterMI function
#' @param grid a vector indicating the grid of values tested for nb.clust. By default 2:5
#' @param graph a boolean indicating if a graphic is plotted
#' @param verbose if TRUE, choosenbclust will print messages on console
#' @param nnodes number of CPU cores for parallel computing. By default, the value used in the call to the clusterMI function
#' @export
#' @importFrom cluster pam agnes
#' @importFrom diceR CSPA
#' @importFrom FactoMineR PCA MCA FAMD
#' @importFrom fpc kmeansCBI nselectboot claraCBI noisemclustCBI hclustCBI
#' @importFrom mix em.mix imp.mix da.mix prelim.mix rngseed
#' @importFrom mclust Mclust
#' @importFrom mice mice complete
#' @importFrom graphics axis
#' @importFrom stats kmeans cutree
#' @importFrom usedist dist_make
#' @seealso \code{\link{imputedata}}
#' @references 
#'  Audigier, V. and Niang, N., Clustering with missing data: which equivalent for Rubin's rules? Advances in Data Analysis and Classification <doi:10.1007/s11634-022-00519-1>, 2022.
#' @examples
#' data(wine)
#' 
#' require(parallel)
#' set.seed(123456)
#' ref <- wine$cult
#' nb.clust <- 3
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- as.matrix(wine.na)
#' wine.na[sample(seq(length(wine.na)), size = ceiling(length(wine.na)/3))] <- NA
#' 
#' # imputation
#' res.imp <- imputedata(data.na=wine.na, nb.clust = nb.clust, m = 5)
#' 
#' # pooling
#' nnodes <- 2 # number of CPU cores for parallel computing
#' res.pool <- clusterMI(res.imp, nnodes = nnodes, instability = FALSE)
#' 
#' # choice of nb.clust
#' \donttest{
#' choosenbclust(res.pool)
#' 
#' }

choosenbclust<-function(output,grid=2:5,graph=TRUE,verbose=TRUE,nnodes=NULL){
  if(verbose){cat("Imputation by ",output$call$output$call$method," for nb.clust between ",min(grid)," and ",max(grid),"...",sep="")}
  if(is.null(nnodes)){nnodes.intern<-output$call$nnodes}else{nnodes.intern<-nnodes}
  
  if(nnodes.intern>1){
    cl <- parallel::makeCluster(nnodes.intern, type = "PSOCK")
    # parallel::clusterEvalQ(cl, library(mice))
    # parallel::clusterEvalQ(cl, library(mclust))
    # parallel::clusterEvalQ(cl, library(mix))
    # parallel::clusterEvalQ(cl, library(fpc))
    parallel::clusterExport(cl, list("grid","output","imputedata","em.mix", "imp.mix", "da.mix", "prelim.mix", "rngseed", "Mclust",
                                     "mice", "complete",
                                     "readData", "createModel", "multipleImp",
                                     "rdirichlet"), envir = environment())
    res.imp<-parallel::parLapply(cl,grid,fun=imputedata,
                    data.na=output$call$output$call$data.na,
                    method=output$call$output$call$method,
                    m=output$call$output$call$m,
                    maxit=output$call$output$call$maxit,
                    Lstart=output$call$output$call$Lstart,
                    L=output$call$output$call$L,verbose=FALSE)
    parallel::stopCluster(cl)
    }else{
  res.imp<-lapply(grid,FUN=imputedata,
                  data.na=output$call$output$call$data.na,
                  method=output$call$output$call$method,
                  m=output$call$output$call$m,
                  maxit=output$call$output$call$maxit,
                  Lstart=output$call$output$call$Lstart,
                  L=output$call$output$call$L,verbose=FALSE)}
  if(verbose){cat(" done!\n")}
  if(verbose){cat("Analysis (by ",output$call$method.clustering,") and pooling (by ",output$call$method.consensus,") ...",sep="")}
  res.rubin<-lapply(res.imp,FUN=clusterMI,method.consensus=output$call$method.consensus,
                    method.clustering=output$call$method.clustering,
                    scaling=output$call$scaling,
                    nb.clust=output$call$nb.clust.intern,
                    Cboot=output$call$Cboot,
                    method.agnes=output$call$method.agnes,
                    modelNames=output$call$modelNames,
                    nstart.kmeans=output$call$nstart.kmeans,
                    m.cmeans=output$call$m.cmeans,
                    nnodes=output$call$nnodes,verbose=FALSE)
 if(verbose){cat(" done!\n")}
  res.out<-sapply(res.rubin,"[[","instability")["Tot",]
  names(res.out)<-grid
  if(graph){
    plot(grid,res.out,xlab="nb clust",ylab="Total instability",type="b",xaxt = "n")
    axis(1, grid, grid)
    }
  res.out<-list(nb.clust=grid[which.min(res.out)],crit=unlist(res.out))
  return(res.out)
}

