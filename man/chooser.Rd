% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chooser.R
\name{chooser}
\alias{chooser}
\title{Kfold cross-validation for specifying threshold r}
\usage{
chooser(
  res.varsel,
  K = 10,
  seed = 12345,
  listvar = NULL,
  grid.r = seq(0, 1, 1/1000),
  graph = TRUE,
  printflag = TRUE,
  nb.clust = NULL,
  nnodes = NULL,
  sizeblock = NULL,
  method.select = NULL,
  B = NULL,
  modelNames = NULL,
  nbvarused = NULL,
  path.outfile = NULL
)
}
\arguments{
\item{res.varsel}{an output from the varselbest function}

\item{K}{an integer given the number of folds}

\item{seed}{a integer}

\item{listvar}{a vector of characters specifiying variables (outcomes) for which cross-validation should be done. By default, all variables that have been considered for varselbest are used.}

\item{grid.r}{a grid for the tuning parameter r}

\item{graph}{a boolean. If TRUE, cross-validation results are printed}

\item{printflag}{a boolean. If TRUE, messages are printed}

\item{nb.clust}{number of clusters. By default, the same as the one used in varselbest}

\item{nnodes}{an integer specifying the number of nodes for parallel computing. By default, the same as  the one used in varselbest}

\item{sizeblock}{number of sampled variables at each iteration. By default, the same as the one used in varselbest}

\item{method.select}{variable selection method used. By default, the same as the one used in varselbest}

\item{B}{number of iterations. By default, the same as the one used in varselbest}

\item{modelNames}{mixture model specification for imputation of subsets. By default, the same as the one used in varselbest}

\item{nbvarused}{a maximal number of selected variables (can be required for a dataset with a large number of variables)}

\item{path.outfile}{a path for message redirection}
}
\value{
A list where each object refers to an outcome variable called in the listvar argument. Each element is composed of three objects
\item{r}{the optimal value for the threshold}
\item{error}{the cross-validation error for each value in \code{grid.r}}
\item{selection}{the subset of selected variables for the optimal threshold}
}
\description{
\code{chooser} returns a list specifying the optimal threshold r for each outcome as well as the associated set of explanatory variables selected, and the cross-validation errror for each value of the grid
}
\details{
\code{varselbest} performs variable selection on random subsets of variables and, then, combines them to recover which explanatory variables are related to the response.
More precisely, the outline of the algorithm are as follows: let consider a random subset of \code{sizeblock} among p variables.
By choosing \code{sizeblock} small, this subset is low dimensional, allowing treatment of missing values by standard imputation method for clustered individuals.
Then, any selection variable scheme can be applied (lasso, stepwise and knockoff are proposed by tuning the \code{method.select} argument).
By resampling \code{B} times, a sample of size \code{sizeblock} among the p variables, we may count how many times, a variable is considered as significantly related to the response and how many times it is not.
We need to define a threshold (\code{r}) to conclude if a given variable is significantly related to the response. \code{chooser} aims at finding the optimal value for the threshold r using Kfold cross-validation.
}
\examples{
data(wine)

require(parallel)
set.seed(123456)
ref <- wine$cult
nb.clust <- 3
wine.na <- wine
wine.na$cult <- NULL
wine.na <- prodna(wine.na)

\donttest{
nnodes <- 2 # parallel::detectCores()
B <- 100 #  Number of iterations
m <- 5 # Number of imputed data sets

# variables selection for incomplete variable "alco"
listvar <- "alco"
res.varsel <- varselbest(data.na = wine.na,
                         nb.clust = nb.clust,
                         listvar = listvar,
                         B = B,
                         nnodes = nnodes)

# frequency of selection
propselect <- res.varsel$proportion[listvar, ]

#predictormatrix with the default threshold value                         
predictmat <- res.varsel$predictormatrix

# r optimal and associated predictor matrix 
res.chooser <- chooser(res.varsel = res.varsel)
thresh <- res.chooser[[listvar]]$r
is.selected <- propselect>=thresh
predictmat[listvar, names(is.selected)] <- as.numeric(is.selected)


# imputation
res.imp.select <- imputedata(data.na = wine.na, method = "FCS-homo",
                     nb.clust = nb.clust, predictmat = predictmat, m = m)
}

}
\references{
Bar-Hen, A. and Audigier, V., An ensemble learning method for variable selection: application to high dimensional data and missing values, Journal of Statistical Computation and Simulation, <doi:10.1080/00949655.2022.2070621>, 2022.
 
 Schafer, J. L. (1997) Analysis of Incomplete Multivariate Data. Chapman & Hall, Chapter 9.
}
