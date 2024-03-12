#' clusterMI: Cluster Analysis with Missing Values by Multiple Imputation
#' 
#' @aliases clusterMI-package
#' 
#' @references
#'  Audigier, V. and Niang, N., Clustering with missing data: which equivalent for Rubin's rules? Advances in Data Analysis and Classification <doi:10.1007/s11634-022-00519-1>, 2022.
#' 
#'  Bar-Hen, A. and Audigier, V., An ensemble learning method for variable selection: application to high dimensional data and missing values, Journal of Statistical Computation and Simulation, <doi:10.1080/00949655.2022.2070621>, 2022.
#' 
#'  Audigier, V., Niang, N., & Resche-Rigon, M. (2021). Clustering with missing data: which imputation model for which cluster analysis method?. arXiv preprint <arXiv:2106.04424>.
#' 
#'  Fang, Y. and Wang, J., Selection of the number of clusters via the bootstrap method. Computational Statistics and Data Analysis, 56, 468-477 <doi:10.1016/j.csda.2011.09.003> 2012.
#' @importFrom graphics abline barplot matplot par plot
#' @importFrom stats coef kmeans lm na.omit step
#' @description
#' \code{clusterMI} is a R package to perform clustering with missing values. For achieving this goal, multiple imputation is used.
#' The package offers various multiple imputation methods dedicated to clustered individuals, as discussed in  Audigier et al. (2021) <arXiv:2106.04424>.
#' In addition, it allows pooling results both in terms of partition and instability, as proposed in Audigier and Niang (2022) <doi:10.1007/s11634-022-00519-1>. Among applications, this instability measure can be used to choose a number of clusters with missing values.
#' @name clusterMI-package
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
#' \donttest{
#' # imputation
#' m <- 5 # number of imputed data sets. Should be larger in practice
#' res.imp <- imputedata(data.na = wine.na, nb.clust = nb.clust, m = m)
#' 
#' # cluster analysis by kmeans and pooling
#' nnodes <- 2 # Number of CPU cores for parallel computing
#' res.pool <- clusterMI(res.imp, nnodes = nnodes)
#' 
#' res.pool$instability
#' table(ref, res.pool$part)
#' 
#' # choice of nb.clust
#' 
#' res.nbclust <- choosenbclust(res.pool)
#' res.nbclust$nb.clust
#' }
#' 
NULL
