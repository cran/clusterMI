#' Kfold cross-validation for specifying threshold r
#'
#' \code{chooser} returns a list specifying the optimal threshold r for each outcome as well as the associated set of explanatory variables selected, and the cross-validation errror for each value of the grid
#' @details \code{varselbest} performs variable selection on random subsets of variables and, then, combines them to recover which explanatory variables are related to the response.
#' More precisely, the outline of the algorithm are as follows: let consider a random subset of \code{sizeblock} among p variables.
#' By choosing \code{sizeblock} small, this subset is low dimensional, allowing treatment of missing values by standard imputation method for clustered individuals.
#' Then, any selection variable scheme can be applied (lasso, stepwise and knockoff are proposed by tuning the \code{method.select} argument).
#' By resampling \code{B} times, a sample of size \code{sizeblock} among the p variables, we may count how many times, a variable is considered as significantly related to the response and how many times it is not.
#' We need to define a threshold (\code{r}) to conclude if a given variable is significantly related to the response. \code{chooser} aims at finding the optimal value for the threshold r using Kfold cross-validation.
#' @return A list where each object refers to an outcome variable called in the listvar argument. Each element is composed of three objects
#' \item{r}{the optimal value for the threshold}
#' \item{error}{the cross-validation error for each value in \code{grid.r}}
#' \item{selection}{the subset of selected variables for the optimal threshold}
#' @export
#' @importFrom mix rngseed
#' @importFrom graphics axis
#' @param res.varsel an output from the varselbest function
#' @param K an integer given the number of folds
#' @param seed a integer
#' @param listvar a vector of characters specifiying variables (outcomes) for which cross-validation should be done. By default, all variables that have been considered for varselbest are used.
#' @param grid.r a grid for the tuning parameter r
#' @param graph a boolean. If TRUE, cross-validation results are printed
#' @param printflag a boolean. If TRUE, messages are printed
#' @param nb.clust number of clusters. By default, the same as the one used in varselbest
#' @param nnodes an integer specifying the number of nodes for parallel computing. By default, the same as  the one used in varselbest
#' @param sizeblock number of sampled variables at each iteration. By default, the same as the one used in varselbest
#' @param method.select variable selection method used. By default, the same as the one used in varselbest
#' @param B number of iterations. By default, the same as the one used in varselbest
#' @param modelNames mixture model specification for imputation of subsets. By default, the same as the one used in varselbest
#' @param nbvarused a maximal number of selected variables (can be required for a dataset with a large number of variables)
#' @param path.outfile a path for message redirection
#' @references 
#'  Bar-Hen, A. and Audigier, V., An ensemble learning method for variable selection: application to high dimensional data and missing values, Journal of Statistical Computation and Simulation, <doi:10.1080/00949655.2022.2070621>, 2022.
#'  
#'  Schafer, J. L. (1997) Analysis of Incomplete Multivariate Data. Chapman & Hall, Chapter 9.
#' 
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
#' \donttest{
#' nnodes <- 2 # parallel::detectCores()
#' B <- 100 #  Number of iterations
#' m <- 5 # Number of imputed data sets
#' 
#' # variables selection for incomplete variable "alco"
#' listvar <- "alco"
#' res.varsel <- varselbest(data.na = wine.na,
#'                          nb.clust = nb.clust,
#'                          listvar = listvar,
#'                          B = B,
#'                          nnodes = nnodes)
#'
#' # frequency of selection
#' propselect <- res.varsel$proportion[listvar, ]
#'
#' #predictormatrix with the default threshold value                         
#' predictmat <- res.varsel$predictormatrix
#' 
#' # r optimal and associated predictor matrix 
#' res.chooser <- chooser(res.varsel = res.varsel)
#' thresh <- res.chooser[[listvar]]$r
#' is.selected <- propselect>=thresh
#' predictmat[listvar, names(is.selected)] <- as.numeric(is.selected)
#'
#' 
#' # imputation
#' res.imp.select <- imputedata(data.na = wine.na, method = "FCS-homo",
#'                      nb.clust = nb.clust, predictmat = predictmat, m = m)
#' }
#' 

chooser<-function(res.varsel,#resultat varselbest
                  K=10,#nb fold
                  seed=12345,
                  listvar=NULL,#nom de la cible
                  grid.r=seq(0,1,1/1000),
                  graph=TRUE,
                  printflag=TRUE,
                  nb.clust=NULL,#vecteur
                  nnodes=NULL,
                  sizeblock=NULL,
                  method.select=NULL,
                  B=NULL,
                  modelNames=NULL,
                  nbvarused=NULL,
                  path.outfile=NULL){
  data.na <- res.varsel$call$data.na
  if(is.null(data.na)){data.na <- res.varsel$call$res.imputedata$call$data.na}
  if(is.null(seed)){stop("a numeric value is expected for seed argument")}
  set.seed(seed)
  rngseed(seed)
  if(is.null(listvar)){listvar<-names(res.varsel$res.varsel)}
  if(is.null(nb.clust)){nb.clust<-sapply(res.varsel$res.varsel[listvar],function(xx){xx$call$nb.clust})}
  if(is.null(sizeblock)){sizeblock <- sapply(res.varsel$res.varsel[listvar],function(xx){xx$call$sizeblock})}
  if(is.null(method.select)){method.select <-sapply(res.varsel$res.varsel[listvar],function(xx){xx$call$methods})}
  if(is.null(B)){B <- sapply(res.varsel$res.varsel[listvar],function(xx){xx$call$B})}
  if(is.null(modelNames)){modelNames <- sapply(res.varsel$res.varsel[listvar],function(xx){xx$call$modelNames})}
  if(is.null(nnodes)){nnodes <- sapply(res.varsel$res.varsel[listvar],function(xx){xx$call$nnodes})}
  if(is.null(nbvarused)){nbvarused<-(ncol(data.na)-1)}
  
  # on effectue le decoupage K fold
  
  indperm<-sample(seq(nrow(data.na)))#permutation des observations
  Kfold<-lapply(seq(K),
                FUN = function(k,indperm,taillek){
                  if(((k+1)*taillek)<=length(indperm)){
                    res<-indperm[seq(from=((k-1)*taillek+1),to=((k)*taillek))]
                  }else{
                    res<-indperm[seq(from=((k-1)*taillek+1),to=length(indperm))]
                  }
                  return(res)
                },
                indperm=indperm,
                taillek=floor(length(indperm)/K)
  )
  if(printflag){cat("Kfold cross-validation using K =",K,"\n") }
  res.out <- mapply(FUN=function(jj,nb.clust,sizeblock,method.select,B,modelNames,nbvarused,grid.r,K,graph,printflag,nnodes,path.outfile,res.varsel){
    if(printflag){cat("Variable ",jj," ")}
    rmin<- -Inf
    erreur<-matrix(NA,K,length(grid.r),dimnames = list(seq(K),grid.r))
    for(k in seq(K)){
      if(printflag){cat(k,"...",sep="")}
      data.train.k<-data.na[-Kfold[[k]],]
      data.test.k<-data.na[Kfold[[k]],]
      res.onefold<-onefold.chooser(data.train.k=data.train.k,
                                   data.test.k=data.test.k,
                                   #res.varsel=res.varsel,#resultat varselbest
                                   jj=jj,#nom de la cible
                                   grid.r=grid.r,
                                   nb.clust=nb.clust,
                                   nnodes=nnodes,
                                   sizeblock=sizeblock,
                                   method.select=method.select,
                                   B=B,
                                   modelNames,K=K,path.outfile = path.outfile,
                                   nbvarused=nbvarused)
      # print(str(res.onefold))
      # save(res.onefold,file=paste0("C:/Users/vince/Desktop/tmp/","res.onefold_",k,".Rdata"))
      res.onefold<-unlist(res.onefold)
      
      rmin<- max(c(as.numeric(names(res.onefold)[1]),rmin))
      
      erreur[k,names(res.onefold)]<- res.onefold
      # save(erreur,file=paste0("C:/Users/vince/Desktop/tmp/","erreur_",k,".Rdata"))
      
      if(graph){
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        par(mfrow=c(1,2))
        plot(grid.r[grid.r>=rmin],
             colSums(erreur[seq(k),grid.r>=rmin,drop=FALSE],na.rm=TRUE),
             xlab="r",
             ylab="error cv",main = paste(c(jj,"error according to r")))
        
        plot(seq(k),
             grid.r[grid.r>=rmin][apply(apply(erreur[,grid.r>=rmin,drop=FALSE],2,cumsum)[seq(k),,drop=FALSE],1,which.min)],
             xlim=c(1,K),
             xlab="Number of test sets",
             ylab="r opt",
             type="b",
             ylim=c(rmin,grid.r[length(grid.r)]),
             main = paste(c(jj,"r opt over folds")),
             axes=FALSE
        )
        axis(1, 1:K, 1:K)
        axis(2)
      }
      
    }
    if(printflag){
      cat("done ! \n")
    }
    erreur<-colSums(erreur[,as.numeric(colnames(erreur))>=rmin],na.rm=TRUE)
    res<-list(r=as.numeric(names(which.min(erreur))),
              error=erreur,
              selection=names(which(res.varsel$proportion[jj,]>as.numeric(names(which.min(erreur))))))
    return(res)
  },
  jj=listvar,nb.clust=nb.clust,sizeblock=sizeblock,method.select=method.select,B=B,modelNames=modelNames,
  MoreArgs = list(grid.r=grid.r,K=K,graph=graph,
                  printflag=printflag,nnodes=nnodes,path.outfile=path.outfile,nbvarused=nbvarused,res.varsel=res.varsel),
  SIMPLIFY = FALSE)
  return(res.out)
}