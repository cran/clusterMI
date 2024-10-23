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
#' data(wine)
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
  if (nnodes.intern > 1) {
    cl <- parallel::makeCluster(nnodes.intern, type = "PSOCK")
    parallel::clusterExport(cl, list("output",
                                     "res.imp",
                                     #fonctions requises par calculintra (fpc)
                                     "kmeansCBI",
                                     "nselectboot", 
                                     "claraCBI",
                                     "noisemclustCBI",
                                     "hclustCBI",
                                     "cmeansCBI.intern",
                                     #fonctions appellees par clusterMI
                                     "fastnmf",
                                     "calculintra.intern",
                                     "cluster.intern",
                                     "clusterMI",
                                     "calcul_inter",
                                     "randindex",
                                     "CSPA",
                                     #fonctions appelees par cluster.intern
                                     "Silhouette.intern",
                                     "cmeans",
                                     "kmeans",
                                     "hclust",
                                     "cutree",
                                     "setTxtProgressBar", 
                                     "txtProgressBar",
                                     "with_seed",
                                     "Mclust",
                                     "hc",
                                     "hcEEE",
                                     "hcVVV", 
                                     "hcVII",
                                     "hcEII",
                                     #factominer
                                     "PCA",
                                     "MCA",
                                     "FAMD"), envir = environment())
    res.rubin <- parallel::parLapply(cl, res.imp, fun = clusterMI,
                                     method.clustering = output$call$method.clustering, 
                                     method.consensus = output$call$method.consensus, 
                                     scaling = output$call$scaling,
                                     nb.clust = output$call$nb.clust.intern, 
                                     Cboot = output$call$Cboot,
                                     method.hclust = output$call$method.hclust,
                                     method.dist =  output$call$method.dist,
                                     modelNames = output$call$modelNames,
                                     modelName.hc=output$call$modelName.hc,
                                     nstart.kmeans = output$call$nstart.kmeans, 
                                     iter.max.kmeans = output$call$iter.max.kmeans,
                                     m.cmeans = output$call$m.cmeans,
                                     nnodes = 1,
                                     instability = TRUE,
                                     verbose = FALSE,
                                     parameter.nmf = output$call$parameter.nmf
                                     )
    parallel::stopCluster(cl)
  }
  else {
    res.rubin <- lapply(res.imp,
                        FUN = clusterMI,
                        method.clustering = output$call$method.clustering, 
                        method.consensus = output$call$method.consensus, 
                        scaling = output$call$scaling,
                        nb.clust = output$call$nb.clust.intern,
                        Cboot = output$call$Cboot,
                        method.hclust = output$call$method.hclust,
                        method.dist =  output$call$method.dist,
                        modelNames = output$call$modelNames,
                        modelName.hc=output$call$modelName.hc,
                        nstart.kmeans = output$call$nstart.kmeans, 
                        iter.max.kmeans = output$call$iter.max.kmeans,
                        m.cmeans = output$call$m.cmeans,
                        nnodes = output$call$nnodes,
                        instability = TRUE,
                        verbose = FALSE,
                        parameter.nmf=output$call$parameter.nmf)
  }
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


