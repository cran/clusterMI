% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmf.R
\name{nmf}
\alias{nmf}
\title{Consensus clustering using non-negative matrix factorization}
\usage{
nmf(
  listpart,
  nb.clust,
  threshold = 10^(-5),
  printflag = TRUE,
  nstart = 100,
  iter.max = 50
)
}
\arguments{
\item{listpart}{a list of partitions}

\item{nb.clust}{an integer specifying the number of clusters}

\item{threshold}{a real specifying when the algorithm is stoped. Default value is 10^(-5)}

\item{printflag}{a boolean. If TRUE, nmf will print messages on console. Default value is TRUE}

\item{nstart}{how many random sets should be chosen for kmeans initalization. Default value is 100}

\item{iter.max}{the maximum number of iterations allowed for kmeans. Default value is 50}
}
\description{
Consensus clustering using non-negative matrix factorization
}
\seealso{
\code{\link[stats]{kmeans}}
}
\keyword{internal}
