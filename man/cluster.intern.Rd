% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.intern.R
\name{cluster.intern}
\alias{cluster.intern}
\title{Apply clustering method after multiple imputation}
\usage{
cluster.intern(
  res.imp,
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
  samples.clara = 500,
  verbose = FALSE
)
}
\arguments{
\item{res.imp}{a list of imputed data sets}

\item{method.clustering}{a single string specifying the clustering algorithm used ("kmeans", "pam", "clara", "hclust" or "mixture","cmeans")}

\item{scaling}{boolean. If TRUE, variables are scaled. Default value is TRUE}

\item{nb.clust}{an integer specifying the number of clusters}

\item{method.hclust}{character string defining the clustering method for hierarchical clustering (required only if method.clustering = "hclust")}

\item{method.dist}{character string defining the method use for computing dissimilarity matrices in hierarchical clustering (required only if method.clustering = "hclust")}

\item{modelNames}{character string indicating the models to be fitted in the EM phase of clustering (required only if method.clustering = "mixture"). By default modelNames = NULL.}

\item{modelName.hc}{A character string indicating the model to be used in model-based agglomerative hierarchical clustering.(required only if method.clustering = "mixture"). By default modelNames.hc = "VVV".}

\item{nstart.kmeans}{how many random sets should be chosen for kmeans initalization. Default value is 100 (required only if method.clustering = "kmeans")}

\item{iter.max.kmeans}{how many iterations should be chosen for kmeans. Default value is 10 (required only if method.clustering = "kmeans")}

\item{m.cmeans}{degree of fuzzification in cmeans clustering. By default m.cmeans = 2}

\item{samples.clara}{number of samples to be drawn from the dataset when performing clustering using clara algorithm. Default value is 500.}

\item{verbose}{logical}
}
\value{
A list with clustering results
}
\description{
From a list of imputed datasets \code{clusterMI} performs cluster analysis on each imputed data set.
}
\details{
Performs cluster analysis (according to the \code{method.clustering} argument). For achieving this goal, the function uses as an input an output from the \code{imputedata} function and applies the cluster analysis method on each imputed data set
 
Step 1 can be tuned by specifying the cluster analysis method used (\code{method.clustering} argument).
If \code{method.clustering = "kmeans"} or \code{"pam"}, then the number of clusters can be specified by tuning the \code{nb.clust} argument. By default, the same number as the one used for imputation is used.
The number of random initializations can also be tuned through the \code{nstart.kmeans} argument.
If \code{method.clustering = "hclust"} (hierarchical clustering), the method used can be specified (see \code{\link[stats]{hclust}}). By default \code{"average"} is used. Furthermore, the number of clusters can be specified, but it can also be automatically chosen if \code{nb.clust} < 0.
If \code{method.clustering = "mixture"} (model-based clustering using gaussian mixture models), the model to be fitted can be tuned by modifying the \code{modelNames} argument (see \code{\link[mclust]{Mclust}}).
If \code{method.clustering = "cmeans"} (clustering using the fuzzy c-means algorithm), then the fuzziness parameter can be modfied by tuning the\code{m.cmeans} argument. By default, \code{m.cmeans = 2}.


Can be performed in parallel by specifying the number of CPU cores (\code{nnodes} argument).
}
\examples{
data(wine)

require(parallel)
set.seed(123456)
ref <- wine$cult
nb.clust <- 3
m <- 5 # number of imputed data sets. Should be larger in practice
wine.na <- wine
wine.na$cult <- NULL
wine.na <- prodna(wine.na)

#imputation
res.imp <- imputedata(data.na = wine.na, nb.clust = nb.clust, m = m)
\donttest{
#analysis by kmeans and pooling
nnodes <- 2 # parallel::detectCores()
res.pool <- clusterMI(res.imp, nnodes = nnodes)

res.pool$instability
table(ref, res.pool$part)
}
}
\seealso{
\code{\link[stats]{hclust}}, \code{\link[mclust]{Mclust}}, \code{\link{imputedata}}, \code{\link[e1071]{cmeans}},\code{\link[stats]{dist}}
}
\keyword{internal}
