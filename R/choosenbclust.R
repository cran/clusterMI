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
#' @importFrom diceR CSPA
#' @importFrom FactoMineR PCA MCA FAMD
#' @importFrom fpc kmeansCBI nselectboot claraCBI noisemclustCBI hclustCBI
#' @importFrom mix em.mix imp.mix da.mix prelim.mix rngseed
#' @importFrom mclust Mclust
#' @importFrom mice mice complete
#' @importFrom graphics axis
#' @importFrom stats kmeans cutree
#' @seealso \code{\link{imputedata}}
#' @references 
#'  Audigier, V. and Niang, N., Clustering with missing data: which equivalent for Rubin's rules? Advances in Data Analysis and Classification <doi:10.1007/s11634-022-00519-1>, 2022.
#' @examples
#' data(wine, package = "clusterMI")
#' 
#' require(parallel)
#' set.seed(123456)
#' ref <- wine$cult
#' nb.clust <- 3
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- prodna(wine.na)
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

choosenbclust<-function (output, grid = 2:5, graph = TRUE, verbose = TRUE, 
                         nnodes = NULL) 
{
  if (verbose) {
    cat("Imputation by ", output$call$output$call$method, 
        " for nb.clust between ", min(grid), " and ", max(grid), 
        "...", sep = "")
  }
  if (is.null(nnodes)) {
    nnodes.intern <- output$call$nnodes
  }
  else {
    nnodes.intern <- nnodes
  }
  if (nnodes.intern > 1) {
    cl <- parallel::makeCluster(nnodes.intern, type = "PSOCK")
    parallel::clusterExport(cl, list("grid", "output",
                                     "imputedata",
                                     #appels imputedata
                                     "myem.mix",
                                     "myimp.mix",
                                     "myimp.cat",
                                     "mclustboot.intern",
                                     "drawW",
                                     "readData",
                                     "createModel", 
                                     "multipleImp",
                                     "rdirichlet",
                                     "rngseedcat",
                                     "em.mix",
                                     "imp.mix",
                                     "da.mix",
                                     "prelim.mix",
                                     "rngseed",
                                    "em.cat",
                                    "imp.cat",
                                    "da.cat",
                                    "prelim.cat",
                                     "Mclust",
                                    "mclustBIC",
                                    "MclustBootstrap",
                                     "mice",
                                     "complete",
                                    "DPMPM_nozeros_imp"), envir = environment())
    res.imp <- parallel::parLapply(cl, grid, fun = imputedata, 
                                   data.na = output$call$output$call$data.na,
                                   method = output$call$output$call$method, 
                                   m = output$call$output$call$m,
                                   maxit = output$call$output$call$maxit, 
                                   Lstart = output$call$output$call$Lstart,
                                   L = output$call$output$call$L,
                                   method.mice=output$call$method.mice,
                                   predictmat = output$call$predictmat,
                                   bootstrap = output$call$bootstrap,
                                   verbose = FALSE)
    parallel::stopCluster(cl)
  }
  else {
    res.imp <- lapply(grid, FUN = imputedata,
                      data.na = output$call$output$call$data.na, 
                      method = output$call$output$call$method,
                      m = output$call$output$call$m, 
                      maxit = output$call$output$call$maxit,
                      Lstart = output$call$output$call$Lstart, 
                      L = output$call$output$call$L,
                      method.mice=output$call$method.mice,
                      predictmat = output$call$predictmat,
                      bootstrap = output$call$bootstrap,
                      verbose = FALSE)
  }
  if (verbose) {
    cat(" done!\n")
  }
  if (verbose) {
    cat("Analysis (by ", output$call$method.clustering, 
        ") and pooling (by ", output$call$method.consensus, 
        ") ...", sep = "")
  }

    res.rubin <- mapply(FUN = function(res.imp,nb.clust,call.clusterMI,nnodes){
      #limitation du temps de calcul de NMF car inutile
      parameter.nmf<- call.clusterMI$parameter.nmf
      parameter.nmf[["threshold"]] <- Inf
      res.out <- clusterMI(res.imp,
                           method.clustering = call.clusterMI$method.clustering,
                           method.consensus = call.clusterMI$method.consensus, 
                           scaling = call.clusterMI$scaling,
                           nb.clust = c("ana"=ifelse(call.clusterMI$nb.clust["ana"]<0,
                                                     call.clusterMI$nb.clust["ana"],
                                                     nb.clust),
                                        "nmf"=nb.clust),
                           Cboot = call.clusterMI$Cboot,
                           method.hclust = call.clusterMI$method.hclust,
                           method.dist =  call.clusterMI$method.dist,
                           modelNames = call.clusterMI$modelNames,
                           modelName.hc=call.clusterMI$modelName.hc,
                           nstart.kmeans = call.clusterMI$nstart.kmeans, 
                           iter.max.kmeans = call.clusterMI$iter.max.kmeans,
                           m.cmeans = call.clusterMI$m.cmeans,
                           nnodes = nnodes,
                           instability = TRUE,
                           verbose = FALSE,
                           parameter.nmf=parameter.nmf
      )
      return(res.out)
    }, res.imp = res.imp , nb.clust = grid, MoreArgs =  list("call.clusterMI" = output$call,"nnodes"=nnodes.intern),SIMPLIFY = FALSE)
  
  if (verbose) {
    cat(" done!\n")
  }
  res.out <- sapply(res.rubin, "[[", "instability")["Tot", 
  ]
  names(res.out) <- grid
  if (graph) {
    plot(grid, res.out, xlab = "nb clust", ylab = "Total instability", 
         type = "b", xaxt = "n")
    axis(1, grid, grid)
  }
  res.out <- list(nb.clust = grid[which.min(res.out)], crit = unlist(res.out))
  return(res.out)
}


