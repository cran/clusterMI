#' Variable selection for specifying conditional imputation models
#' 
#' @description
#' \code{varselbest} performs variable selection from an incomplete dataset (see Bar-Hen and Audigier (2022) <doi:10.1080/00949655.2022.2070621>) in order to specify the imputation models to use for FCS imputation methods
#' @details
#' \code{varselbest} performs variable selection on random subsets of variables and, then, combines them to recover which explanatory variables are related to the response.
#' More precisely, the outline of the algorithm are as follows: let consider a random subset of \code{sizeblock} among p variables.
#' By choosing \code{sizeblock} small, this subset is low dimensional, allowing treatment of missing values by standard imputation method for clustered individuals.
#' Then, any selection variable scheme can be applied (lasso, stepwise and knockoff are proposed by tuning the \code{method.select} argument).
#' By resampling \code{B} times, a sample of size \code{sizeblock} among the p variables, we may count how many times, a variable is considered as significantly related to the response and how many times it is not.
#' We need to define a threshold (\code{r}) to conclude if a given variable is significantly related to the response.
#' @return a list of four objects
#' \item{predictormatrix}{a numeric matrix containing 0 and 1 specifying on each line the set of predictors to be used for each target column of the incomplete dataset.}
#' \item{res.varsel}{a list given details on the variable selection procedure (only required for checking convergence by the \code{chooseB} function)}
#' \item{proportion}{a numeric matrix of proportion indicating on each line the variable importance of each predictor}
#' \item{call}{the matching call}
#' @param data.na a dataframe with only numeric variables
#' @param res.imputedata an output from \code{\link{imputedata}}
#' @param listvar a character vector indicating for which subset of incomplete variables variable selection must be performed. By default all column names.
#' @param nb.clust the number of clusters used for imputation
#' @param nnodes number of CPU cores for parallel computing. By default, nnodes = 1
#' @param sizeblock an integer indicating the number of variables sampled at each iteration
#' @param method.select a single string indicating the variable selection method applied on each subset of variables
#' @param B number of iterations, by default B = 200
#' @param r a numerical vector (or a single real number) indicating the threshold used for each variable in listvar. Each value of r should be between 0 and 1. See details.
#' @param graph a boolean. If TRUE two graphics are plotted per variable in \code{listvar}: a graphic reporting the variable importance measure of each explanatory variable and a graphic reporting the influence of the number iterations (B) on the importance measures
#' @param printflag a boolean. If TRUE, a message is printed at each iteration. Use printflag = FALSE for silent selection.
#' @param path.outfile a vector of strings indicating the path for redirection of print messages. Default value is NULL, meaning that silent imputation is performed. Otherwise, print messages are saved in the files path.outfile/output.txt. One file per node is generated.
#' @param mar a numerical vector of the form c(bottom, left, top, right). Only used if graph = TRUE
#' @param cex.names expansion factor for axis names (bar labels) (only used if graph = TRUE)
#' @param modelNames a vector of character strings indicating the models to be fitted in the EM phase of clustering
#' @seealso \code{\link[mice]{mice}},\code{\link{clusterMI}}, \code{\link{imputedata}},\code{\link[knockoff]{knockoff}},\code{\link[glmnet]{glmnet}},\code{\link[mix]{imp.mix}}
#' @export
#' @importFrom mice mice
#' @importFrom graphics title
#' @references 
#'  Bar-Hen, A. and Audigier, V., An ensemble learning method for variable selection: application to high dimensional data and missing values, Journal of Statistical Computation and Simulation, <doi:10.1080/00949655.2022.2070621>, 2022.
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
#' nnodes <- 2 # parallel::detectCores()
#' B <- 150 #  Number of iterations
#' m <- 5 # Number of imputed data sets
#' 
#' # variable selection
#' res.varsel <- varselbest(data.na = wine.na,
#'                          nb.clust = nb.clust,
#'                          listvar = c("alco","malic"),
#'                          B = B,
#'                          nnodes = nnodes)
#' predictmat <- res.varsel$predictormatrix
#' 
#' # imputation
#' res.imp.select <- imputedata(data.na = wine.na, method = "FCS-homo",
#'                      nb.clust = nb.clust, predictmat = predictmat, m = m)
#' }
#' 

varselbest<-function(data.na=NULL,res.imputedata=NULL,listvar=NULL,nb.clust=NULL,nnodes=1,sizeblock = 5,
                     method.select="knockoff",B=200,r=.3,graph=TRUE, printflag=TRUE, path.outfile=NULL,
                     mar=c(2, 4, 2, .5) + 0.1,cex.names=.7, modelNames=NULL
                     ){
  Call <- list(data.na=data.na,
               res.imputedata=res.imputedata,
               listvar=listvar,
               nb.clust=nb.clust,
               nnodes=nnodes,
               sizeblock = sizeblock,
               method.select=method.select,
               B=B,
               r=r,
               graph=graph,
               printflag=printflag,
               path.outfile=path.outfile,
               mar=mar,
               cex.names=cex.names,
               modelNames=modelNames)
  
  if(is.null(data.na)&is.null(res.imputedata)){stop("data.na or res.imputedata should be given")}
  if(is.null(nb.clust)&is.null(res.imputedata)){stop("nb.clust or res.imputedata should be given")}
if(is.null(data.na)){
data.na.intern<-res.imputedata$call$data.na
}else{data.na.intern<-data.na
if(sizeblock>=ncol(data.na)){stop("sizeblock should be less than the number of variables")}}
  if(any(sapply(data.na.intern,is.factor))|any(sapply(data.na.intern,is.ordered))){
    stop("data must be continuous")
  }
  
  
if(is.null(nb.clust)){nb.clust.intern<-res.imputedata$call$nb.clust}else{nb.clust.intern<-nb.clust}
# if(is.null(nnodes)){nnodes.intern<-res.imputedata$call$nnodes}else{nnodes.intern<-nnodes}
nnodes.intern<-nnodes
if(is.null(listvar)){listvar.intern<-colnames(data.na.intern)}else{listvar.intern<-listvar}
if(length(r)==length(listvar.intern)){r.intern<-r}else if(length(r)==1){r.intern<-rep(r,length(listvar.intern))}else{stop("non-conforming size for r")}
names(r.intern)<-listvar.intern
if(graph){
  Mfrow_old <- par()$mfrow
  Mfrow <- c(min(length(listvar.intern), 4), 1 + (length(listvar.intern) - 1)%/%4)
  par(mfrow = Mfrow,mar=mar)
}

propmatrix<-predictormatrix<-mice(data = data.na.intern,maxit=0,seed=1234)$predict
res.varsel<-list()

for (varindex in listvar.intern){
  if (printflag){cat(varindex,"...",sep="")}
  oldname<-colnames(data.na.intern)
  names(oldname)<-oldname
  
  oldname[varindex]<-"Y"
  oldname[-which(oldname=="Y")]<-paste0("V",seq(ncol(data.na.intern)-1))
  X<-as.matrix(data.na.intern[,-which(varindex==names(oldname))])
  Y<-matrix(data.na.intern[,names(which(oldname=="Y"))],ncol=1)
  colnames(X)<-oldname[-which(oldname=="Y")]
  colnames(Y)<-"Y"
  ind.full<-which(!(is.na(Y)))
  
  res.varsel[[varindex]]<-varselbest.intern(nnodes=nnodes.intern,
                                            X =X[ind.full,],
                                            Y = Y[ind.full,drop=FALSE],
                                            B=B,sizeblock = sizeblock,methods=method.select,
                                            path.outfile=path.outfile,nb.clust = nb.clust.intern,
                                            printflag=FALSE,r=r.intern[varindex],modelNames=modelNames
  )
  tmp.bar<-res.varsel[[varindex]]$res$proportion
  names(tmp.bar)<-names(oldname)[-which(names(oldname)==varindex)]
  res.varsel[[varindex]]$res$proportion<-tmp.bar
 
  if(graph){
      barplot(sort(tmp.bar,decreasing=TRUE), ylab="proportion",
              main=varindex,ylim=c(0,1),las=2,cex.names = cex.names)
      abline(h=r.intern[varindex],col=2,lty=2)
  }
  
  #rename
  var.in<-names(oldname)[which(oldname%in%res.varsel[[varindex]]$res$selection)]
  res.varsel[[varindex]]$res$selection<-var.in
  
  names(res.varsel[[varindex]]$res$garde)<-names(oldname)[which(oldname%in%names(res.varsel[[varindex]]$res$garde))]
  names(res.varsel[[varindex]]$res$effectif)<-names(oldname)[which(oldname%in%names(res.varsel[[varindex]]$res$effectif))]
  names(res.varsel[[varindex]]$res$failure)<-names(oldname)[which(oldname%in%names(res.varsel[[varindex]]$res$failure))]
  
  # var.out<-names(oldname)[which(oldname%in%names(which(res.varsel[[varindex]]$res$proportion<r.intern[varindex])))]
  var.out<-names(which(res.varsel[[varindex]]$res$proportion<r.intern[varindex]))
  predictormatrix[varindex,var.out]<-0
  propmatrix[varindex,names(tmp.bar)]<-tmp.bar
}


res.out<-list(predictormatrix=predictormatrix,
              res.varsel=res.varsel,
              proportion=propmatrix,
              call=Call)
if(graph){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = Mfrow,mar=mar)
  tmp<-mapply(res.varsel,listvar.intern,FUN = function(xx,varindex){chooseB.intern(xx);title(varindex)})
  # par(mfrow = Mfrow_old)
}
if (printflag){cat("done!\n")}
return(res.out)
}