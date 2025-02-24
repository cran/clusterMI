#' Multiple imputation methods for cluster analysis
#'
#' \code{imputedata} returns a list of imputed datasets by using imputation methods dedicated to individuals clustered in (unknown) groups
#'
#' @details The \code{imputedata} offers various multiple imputation methods dedicated to clustered individuals.
#' In particular, two fully conditional imputation methods are proposed (\code{FCS-homo} and \code{FCS-hetero}) which essentially differ by the assumption about the covariance in each cluster (constant or not respectively).
#' The imputation requires a pre-specified number of clusters (\code{nb.clust}). See \code{\link{choosenbclust}} if this number is unknown.
#' The \code{imputedata} function alternates clustering and imputation given the partition of individuals.
#' When the clustering is performed, the function calls the \code{\link[mice]{mice}} function from the \code{mice} R package to perform imputation.
#' The \code{mice} package proposes various methods for imputation which can be specified by tuning the \code{method.mice} argument.
#' Note that two other joint modelling methods are also available: \code{JM-GL} from the R package \code{mix} and \code{JM-DP} from the R package \code{DPImputeCont} \url{https://github.com/hang-j-kim/DPImputeCont}
#' @return a list of 3 objets
#'  \item{res.imp}{ a list with the several imputed datasets}
#'  \item{res.conv}{ for FCS methods, an array given the between (and within) inertia of each imputed variable at each iteration and for each imputed dataset. For JM methods, a matrix given the between inertia for each variable and each imputed dataset.}
#'  \item{call}{ the matching call}
#' 
#' @param data.na an incomplete dataframe
#' @param method a single string specifying the imputation method used among "FCS-homo","FCS-hetero","JM-DP","JM-GL". By default method = "JM-GL". See the details section
#' @param nb.clust number of clusters
#' @param m number of imputed datasets. By default, m = 20.
#' @param maxit number of iterations for FCS methods (only used for method = FCS-homo or method = FCS-hetero)
#' @param Lstart number of iterations for the burn-in period (only used if method ="JM-DP" or "JM-GL")
#' @param L number of skipped iterations to keep one imputed data set after the burn-in period (only used if method ="JM-DP" or "JM-GL")
#' @param method.mice a vector of strings (or a single string) giving the imputation method for each variable (only used for method = FCS-homo or method = FCS-hetero). Default value is "pmm" (predictive mean matching) for FCS-homo and "mice.impute.2l.jomo" for FCS-hetero
#' @param predictmat predictor matrix used for FCS imputation (only used for method = FCS-homo or method = FCS-hetero)
#' @param verbose a boolean. If TRUE, a message is printed at each iteration. Use verbose = FALSE for silent imputation
#' @param seed a positive integer initializing the random generator
#' @param bootstrap a boolean. Use bootstrap = TRUE for proper imputation with FCS methods (Mclust sometimes fails with multiple points)
#' @seealso \code{\link[mice]{mice}} \code{\link{choosenbclust}} \code{\link{choosemaxit}} \code{\link{varselbest}} \code{\link[mix]{imp.mix}}
#' @export
#' @importFrom mix em.mix imp.mix da.mix prelim.mix rngseed
#' @importFrom cat em.cat imp.cat da.cat prelim.cat
#' @rawNamespace importFrom(cat, rngseedcat=rngseed)
#' @importFrom FactoMineR MCA
#' @importFrom mclust Mclust mclustBIC MclustBootstrap
#' @importFrom mice mice complete
#' @importFrom NPBayesImputeCat DPMPM_nozeros_imp
#' @references
#' Kim, H. J., Reiter, J. P., Wang, Q., Cox, L. H. and Karr, A. F. (2014), Multiple imputation of missing or faulty values under linear constraints, Journal of Business and Economics Statistics, 32, 375-386 <doi:10.1080/07350015.2014.885435>
#' 
#' Schafer, J. L. (1997) Analysis of Incomplete Multivariate Data. Chapman & Hall, Chapter 9.
#' 
#' Audigier, V., Niang, N., & Resche-Rigon, M. (2021). Clustering with missing data: which imputation model for which cluster analysis method?. arXiv preprint <arXiv:2106.04424>.
#' @examples
#' data(wine, package = "clusterMI")
#' set.seed(123456)
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- prodna(wine.na)
#' nb.clust <- 3 # number of clusters
#' m <- 20 # number of imputed data sets
#' res.imp <- imputedata(data.na = wine.na, nb.clust = nb.clust, m = m)
#' lapply(res.imp$res.imp, summary)

imputedata<-function(data.na,
                     method="JM-GL",
                    nb.clust=NULL,
                    m=20,
                    maxit=50,
                    Lstart=100,
                    L=20,
                    method.mice=NULL,
                    predictmat=NULL,
                    verbose=TRUE,
                    seed=1234,
                    bootstrap = FALSE){
  if(!is.data.frame(data.na)){stop("data.na need to be a dataframe.")}
  if(sum(is.na(data.na))==0){stop("The dataset must be incomplete")}
  if(is.null(nb.clust)){stop("The number of clusters (nb.clust) must be specified")}
  if(is.null(colnames(data.na))){stop("Colnames of data.na are missing")}
  data.type<-"continuous"
  if(any(sapply(data.na,is.factor))|any(sapply(data.na,is.ordered))){
    data.type<-"categorical"
    if(any(sapply(data.na,is.numeric))|any(sapply(data.na,is.integer))){
      data.type<-"mixed"
    }
  }
  # check number of categories for each categorical variable
  if(data.type != "continuous"){
    tmp<-(sapply(data.na[sapply(data.na,is.factor)],nlevels)==1)
    if(any(tmp)){warning(cat("At least one factor has only one level and cannot be considered for clustering. Variable(s)",names(which(tmp)),"should be removed.\n"))}
  }
  
  #check colnames (mice ne supporte pas n'importe quel nom)
  checkcolnames<-try(sapply(colnames(data.na),str2lang),silent=TRUE)
  if(inherits(checkcolnames,"try-error")&method%in%c("FCS-homo","FCS-hetero")){
    stop("Colnames of data.na are not supported by mice. Please rename colnames with character string compatible with str2lang")
  }
  
  #check colnames for predicmat
  if(!is.null(predictmat)){
    if(!(identical(colnames(predictmat),colnames(data.na)))){
      stop("The column names for predictmat and data.na do not match")}
  }
  
  
  Call<-list(data.na=data.na,
             method=method,
  nb.clust=nb.clust,
  m=m,
  maxit=maxit,
  Lstart=Lstart,
  L=L,
  method.mice=method.mice,
  predictmat=predictmat,
  verbose=TRUE,
  seed=1234,
  data.type=data.type,
  bootstrap=bootstrap)
  
  if(is.null(method.mice)){
    method.mice.homo<-NULL;method.mice.hetero<-"2l.jomo"
  }else{
    if(length(method.mice)==1){
      method.mice.hetero<- method.mice.homo<-method.mice
    }else{
      if("FCS-homo"%in%method){
        method.mice.homo<-c("",method.mice)
      }else if("FCS-hetero"%in%method){
        method.mice.hetero<-c("",method.mice)}
    }
  }

  res.conv<-NULL
  
  if(verbose){cat(method,"\n")}
  
  if(("FCS-homo"== method)){
    if(data.type=="continuous"){
      data.na.schaf<-data.na
      prior.dir<-try(myem.mix(dataset = data.na,
                              K = nb.clust,silent = TRUE)$theta$pi,silent=TRUE)#ML
      quali.index<-NULL
      res.conv<-array(NA,dim=c(m,maxit,ncol(data.na),2),dimnames=list(seq.int(m),seq.int(maxit),colnames(data.na),c("var_inter","var_intra")))
      
    }else if(data.type=="categorical"){
      data.na.schaf<-sapply(data.na,as.numeric)
      prior.dir<-0.5;class(prior.dir)<-"try-error"
      res.conv<-NULL
    }else if(data.type=="mixed"){
      quali.index<-which(sapply(data.na,is.factor)|(sapply(data.na,is.ordered)))
      quanti.index<-which(!(sapply(data.na,is.factor)|(sapply(data.na,is.ordered))))
      
      data.na.schaf<-data.na
      data.na.schaf[,quali.index]<-sapply(data.na.schaf[,quali.index],as.numeric)
      data.na.schaf<-data.na.schaf[,c(quali.index,quanti.index)]
      prior.dir<-0.5;class(prior.dir)<-"try-error"
      res.conv<-array(NA,dim=c(m,maxit,ncol(data.na)-length(quali.index),2),
                      dimnames=list(seq.int(m),seq.int(maxit),colnames(data.na)[quanti.index],c("var_inter","var_intra")))
      
    }
    if(inherits(prior.dir,"try-error")){
      prior.dir<-0.5
    }
    res.imp<-list()
    for (tab in seq.int(m)) {
      seed.tmp <- NULL
      if (tab == 1) {
        seed.tmp <- seed
      }
      if (verbose) {
        cat("\n m =", tab)
      }
      flag <- FALSE
      if (data.type %in% c("continuous", "mixed")) {
        data.init <- try(as.data.frame(myimp.mix(data.na.schaf, 
                                                 K = nb.clust, silent = TRUE, seed = seed.tmp, 
                                                 nbquali = length(quali.index))),silent=TRUE)
        if ((inherits(data.init,"try-error")) & (tab == 1)) {
          mix::rngseed(seed)
        }
      }
      else if (data.type == "categorical") {
        data.init <- try(as.data.frame(myimp.cat(data.na.schaf, 
                                                 K = nb.clust, silent = TRUE, seed = seed.tmp)),silent = TRUE)
        if ((inherits(data.init,"try-error")) & (tab == 1)) {
          rngseedcat(seed)
        }
      }
      if (!inherits(data.init,"try-error")) {
        data.init <- as.data.frame(data.init)
        for (variable in names(quali.index)) {
          data.init[[variable]] <- factor(data.init[[variable]], 
                                          labels = levels(data.na[[variable]]))
        }
        if (is.null(predictmat)) {
          if (min(table(data.init$class)) < ncol(data.init)) {
            flag <- TRUE
          }
        }
        else if (!is.null(predictmat)){
          if (min(table(data.init$class)) <= min(rowSums(predictmat))) {
            flag <- TRUE
          }
        }
      }
      if ((!inherits(data.init,"try-error")) & (!flag)) {
        flag <- (length(table(data.init$class)) != nb.clust)
      }
      if (inherits(data.init,"try-error") | flag){
        if(verbose){cat("random initialization\n")}
        data.init <- as.data.frame(data.na)
        data.init <- complete(mice(data = data.init, 
                                   maxit = 0))
        clustering.tmp <- sample(seq.int(nb.clust), size = nrow(data.na), 
                                 replace = TRUE)
        data.init.class <- cbind.data.frame(class = as.factor(clustering.tmp), 
                                            data.init)
      }else{
        clustering.tmp<-data.init$class
        data.init.class<-cbind.data.frame(class=as.factor(clustering.tmp),
                                          data.init[,-which(colnames(data.init)=="class")])
      }
      
      data.na.class<-cbind.data.frame(class=as.factor(clustering.tmp),data.na)
      if(verbose){cat(" nbiter = ")}
      for(nbiter in seq.int(maxit)){
        if(verbose){cat(nbiter,"...",sep="")}
        #besoin de mettre a jour le clustering a part, car completement manquant
        if(is.null(predictmat)){
          # on vérifie l'effectif des classes en fonction du nombre de variables
          if (min(table(data.na.class$class)) < ncol(data.na.class)) {
            warning("The number of individuals per cluster is too low to build the regression model. The imputation cycle is stopped.")
            break()
          }
          
          mice.tmp <- mice(data = data.na.class, m = 1, 
                           method = method.mice.homo, maxit = 1, printFlag = FALSE, 
                           data.init = data.init.class)
          if (!is.null((mice.tmp$loggedEvents))) {
            warning("mice returns a potential issue\n", mice.tmp$loggedEvents) 
          }
          
          data.init.class<-complete(mice.tmp)
          
        }else if(!is.null(predictmat)){
          # on vérifie l'effectif des classes en fonction du nombre de variables
          if (min(table(data.na.class$class)) <= min(rowSums(predictmat))) {
            warning("The number of individuals per cluster is too low to build the regression model. The imputation cycle is stopped.")
            break()
          }
          predictormatrix<-rbind(c(0,rep(0,nrow(predictmat))),cbind(rep(1,nrow(predictmat)),predictmat))
          colnames(predictormatrix)<-c("class",colnames(predictmat))
          rownames(predictormatrix)<-c("class",rownames(predictmat))
          data.init.class<-complete(mice(data = data.na.class,
                                         m = 1,
                                         method = method.mice.homo,
                                         maxit=1,
                                         printFlag=FALSE,
                                         data.init = data.init.class,
                                         predictorMatrix = predictormatrix))
        }
        
        varnotimputed<-names(which(apply(is.na(data.init.class),2,any)))
        if(length(varnotimputed)>0){warning(c("The following variables are not imputed: ",varnotimputed,". This can be due to colinearity between variables."))}
        if(data.type=="continuous"){
          data.init.class.intern<-data.init.class
          if(length(varnotimputed>0)){data.init.class.intern<-data.init.class[,-which(colnames(data.init.class)%in%varnotimputed)]}
          
        }else if(data.type=="categorical"){
          
          data.init.class.intern<-cbind.data.frame(class=data.init.class[,which(colnames(data.init.class)=="class")],
                                                   MCA(
                                                     data.init.class[,-which(colnames(data.init.class)%in%c("class",varnotimputed))],
                                                     graph=FALSE,ncp=Inf)$ind$coord)
          
        }else if(data.type=="mixed"){
          if(length(intersect(names(quali.index),varnotimputed))==length(quali.index)){
            #il n'y a alors plus de variables quali
            data.init.class.intern<-cbind.data.frame(class=data.init.class[,which(colnames(data.init.class)=="class")],
                                                     data.init.class[,-which(colnames(data.init.class)%in%c("class",varnotimputed))])
            
          }else{
            data.init.class.intern<-cbind.data.frame(class=data.init.class[,which(colnames(data.init.class)=="class")],
                                                     FAMD(
                                                       data.init.class[,-which(colnames(data.init.class)%in%c("class",varnotimputed))],
                                                       graph=FALSE,ncp=Inf)$ind$coord)
          }
        }
        if(!bootstrap){
          # tirage de W
          res.myem <- myem.mix(dataset = data.init.class.intern[, 
                                                                -which(colnames(data.init.class.intern) == 
                                                                         "class"),drop = FALSE],
                               K = nb.clust,
                               silent = TRUE)
          res.myimp <- myimp.mix(data.init.class.intern[, 
                                                        -which(colnames(data.init.class.intern) == 
                                                                 "class"), drop=FALSE], K = nb.clust, seed = NULL, 
                                 silent = TRUE)
          W<-res.myimp[,"class"]
          # tirage de theta
          res.myem.new<-res.myem
          res.myem.new$s$w<-W
          newtheta<-da.mix(res.myem.new$s,res.myem.new$theta,steps=1,prior=prior.dir)
          res.myem.new$theta$pi<-newtheta$pi
          res.myem.new$s$w<-res.myem$s$w
          
          #imputation
          
          W<-imp.mix(res.myem.new$s,theta = res.myem.new$theta)[,"class"]
          
        }else if(bootstrap){
          if(ncol(data.init.class.intern[,-which(colnames(data.init.class.intern)=="class"),drop=FALSE])==1){
            #une seule variable pour le clustering
            modelNames.tmp<-"E"
          }else{
            modelNames.tmp<- "EEE"
          }
          res.mclust<-mclustboot.intern(data.init.class.intern[,-which(colnames(data.init.class.intern)=="class"),drop=FALSE],
                                        G=nb.clust,modelNames = modelNames.tmp, verbose =FALSE)
          res.mclust$data<-as.matrix(data.init.class.intern[,-which(colnames(data.init.class.intern)=="class"),drop=FALSE])
          
          W<-drawW(res.mclust = res.mclust,method = "LDA")
          #gestion du cas ou une classe est vide
          if(length(unique(W))<nb.clust){
            tableW<-rep(0,nb.clust)
            names(tableW)<-as.character(seq.int(nb.clust))
            tableW[as.character(unique(W))]<-table(W)
          }else{
            tableW<-table(W)
          }
          
          
          theta<-as.vector(rdirichlet(1,alpha = prior.dir+tableW))
          res.mclust$parameters$pro<-theta
          W<-drawW(res.mclust = res.mclust,method = "LDA")
        }
        clustering.tmp<-W
        data.na.class<-cbind.data.frame(class=as.factor(as.character(W)),data.na)
        data.init.class$class<-as.factor(as.character(W))
        if(data.type%in%c("mixed","continuous")){
          res.conv[tab, nbiter, , ] <- t(apply(
            data.init.class[, -which(colnames(data.init.class) %in% c(names(quali.index), "class")), drop = FALSE],
            2,
            FUN=function(xx,ref){
              res.by<-do.call(rbind,by(xx,INDICES=ref,
                                       FUN=function(yy){
                                         meangp<-mean(yy)
                                         fregp<-length(yy)
                                         res.out<-c(mean=meangp,freq=fregp)
                                         return(res.out)
                                       }))
              inter<-sum(((res.by[,"mean"]-mean(xx))^2)*(res.by[,"freq"]))/sum(res.by[,"freq"])
              intra<-mean((xx-mean(xx))^2)-inter
              res.out<-c("inter"=inter,"intra"=intra)
              return(res.out)
            },ref=as.factor(data.init.class$class)))
        }
      }
      
      res.imp[[tab]]<-data.init.class
      res.imp[[tab]]$class<-NULL
    }
    if(verbose){cat("done!\n")}
  }
 
  if(("FCS-hetero"== method)){
    if(data.type!="continuous"){stop("FCS-hetero is not available for non-continuous data. Use FCS-homo (or a JM method).")}
    quali.index<-NULL
       prior.dir<-try(myem.mix(dataset = data.na,
                        K = nb.clust,silent = TRUE)$theta$pi,silent=TRUE)
    if("try-error"%in%class(prior.dir)){prior.dir<-0.5}
    res.conv<-array(NA,dim=c(m,maxit,ncol(data.na),2),dimnames=list(seq.int(m),seq.int(maxit),colnames(data.na),c("var_inter","var_intra")))
    res.imp<-list()
    # if(checkconv){par(mfrow=c(m,1))}
    for(tab in seq.int(m)){
      if(verbose){cat("\n m =",tab)}
      if(tab==1){data.init<-try(as.data.frame(myimp.mix(data.na,K = nb.clust,silent = TRUE, seed = seed)),silent=TRUE)
      }else{
        data.init<-try(as.data.frame(myimp.mix(data.na,K = nb.clust,seed=NULL,silent = TRUE)),silent=TRUE)
      }
      if("try-error"%in%class(data.init)){
        if(verbose){cat("random initialization\n")}
        data.init<-as.data.frame(data.na)
        data.init<-cbind.data.frame(class=as.factor(sample(seq.int(nb.clust),size=nrow(data.na),replace=TRUE)),complete(mice(data = data.init,maxit=0)))
      }
      clustering.tmp<-as.integer(data.init$class)
      data.init.tmp<-cbind.data.frame(class=clustering.tmp,data.init[,-which(colnames(data.init)=="class")])
      
      
      don.na.class<-cbind.data.frame(class=clustering.tmp,data.na)
      ind.clust<-1
      
      if(is.null(predictmat)){
        predictor.matrix<-mice(cbind.data.frame(class=clustering.tmp,
                                                data.na),m=1,maxit=0)$pred
        predictor.matrix[ind.clust,ind.clust]<-0
        predictor.matrix[-ind.clust,ind.clust]<- -2
        predictor.matrix[predictor.matrix==1]<-2
      }else{
        predictor.matrix<-rbind(c(0,rep(0,nrow(predictmat))),cbind(rep(1,nrow(predictmat)),predictmat))
        colnames(predictor.matrix)<-c("class",colnames(predictmat))
        rownames(predictor.matrix)<-c("class",rownames(predictmat))
        predictor.matrix[ind.clust,ind.clust]<-0
        predictor.matrix[-ind.clust,ind.clust]<- -2
        predictor.matrix[predictor.matrix==1]<-2
      }
      method.mice.tmp<-rep(method.mice.hetero,ncol(data.init));names(method.mice.tmp)<-colnames(data.init);method.mice.tmp[method=="class"]<-""
      if(verbose){cat(" nbiter = ")}
      for(nbiter in seq.int(maxit)){
        if(verbose){cat(nbiter,"...",sep="")}
        #besoin de mettre a jour le clustering ? part, car completement manquant
        data.init.tmp<-complete(mice(data = don.na.class,
                                     m = 1,
                                     method = method.mice.tmp,
                                     maxit=1,
                                     printFlag=FALSE,
                                     data.init = data.init.tmp,
                                     predictorMatrix = predictor.matrix))
        
        varnotimputed<-names(which(apply(is.na(data.init.tmp),2,any)))
        if(length(varnotimputed)>0){warning(c("The following variables are not imputed: ",varnotimputed,". This can be due to colinearity between variables."))}
        
        # tirage de W
        if(!bootstrap){
          res.init.hc <- hc(data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))],
                            modelName = "VVV",
                            use = "STD")
          res.mclust<-try(Mclust(data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))],
                             G=nb.clust,
                             initialization = list(hcPairs = res.init.hc),
                             modelNames = "VVV",
                             verbose = FALSE),silent=TRUE)
          if((is.null(res.mclust))|("try-error"%in%class(res.mclust))){
            res.mclust<-Mclust(data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))],
                               G=nb.clust,
                               initialization = list(hcPairs = res.init.hc),
                               modelNames = "EEE",verbose = FALSE)
            warning("Model too complex, EEE is fitted")
          }
        }else{
          res.mclust<-try(mclustboot.intern(data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))],
                                        G=nb.clust,modelNames = "VVV",verbose = FALSE),silent=TRUE)
          if((is.null(res.mclust))|("try-error"%in%class(res.mclust))){
            res.mclust<-mclustboot.intern(data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))],
                                          G=nb.clust,modelNames = "EEE",verbose = FALSE)
            warning("Model too complex, EEE is fitted")
            if(is.null(res.mclust)){stop("Model EEE cannot be fitted")}
          }
          res.mclust$data<-as.matrix(data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))])
        }
        Wstar<-drawW(res.mclust, method = "QDA")
        
        #gestion du cas ou une classe est vide
        if(length(unique(Wstar))<nb.clust){
          tableWstar<-rep(0,nb.clust)
          names(tableWstar)<-as.character(seq.int(nb.clust))
          tableWstar[as.character(unique(Wstar))]<-table(Wstar)
        }else{
          tableWstar<-table(Wstar)
        }
        # tirage de theta
        theta<-rdirichlet(n = 1,alpha = prior.dir+tableWstar)
        #imputation
        res.mclust$parameters$pro<-as.vector(theta)
        W<-drawW(res.mclust, method = "QDA")
        clustering.tmp<-as.factor(W)
        don.na.class<-cbind.data.frame(class=as.integer(W),data.na)
        data.init.tmp$class<-clustering.tmp
        
        res.conv[tab, nbiter, , ] <- t(apply(
          data.init.tmp[, -which(colnames(data.init.tmp) %in% c(names(quali.index), "class")), drop = FALSE],
          2,
          FUN=function(xx,ref){
            res.by<-do.call(rbind,by(xx,INDICES=ref,
                                     FUN=function(yy){
                                       meangp<-mean(yy)
                                       fregp<-length(yy)
                                       res.out<-c(mean=meangp,freq=fregp)
                                       return(res.out)
                                     }))
            inter<-sum(((res.by[,"mean"]-mean(xx))^2)*(res.by[,"freq"]))/sum(res.by[,"freq"])
            intra<-mean((xx-mean(xx))^2)-inter
            res.out<-c("inter"=inter,"intra"=intra)
            return(res.out)
          },ref=as.factor(data.init.tmp$class)))
        
        
      }
      res.imp[[tab]]<-data.init.tmp
      res.imp[[tab]]$class<-NULL
    }
    if(verbose){cat("done!\n")}
  }
 
  if("JM-GL"== method){
    res.imp<-list()
    rngseed(seed)
    temp.part<-matrix(NA,nrow=m,ncol=nrow(data.na))
    if(data.type%in%c("continuous","mixed")){
     
      data.na.intern<-data.na
        quali.index<-which(sapply(data.na,is.factor)|(sapply(data.na,is.ordered)))
        quanti.index<-which(!(sapply(data.na,is.factor)|(sapply(data.na,is.ordered))))
        if(data.type=="mixed"){
          quanti.index<-which(!(sapply(data.na,is.factor)|(sapply(data.na,is.ordered))))
          data.na.intern[,quali.index]<-sapply(data.na.intern[,quali.index],as.numeric)
          data.na.intern<-data.na.intern[,c(quali.index,quanti.index)]
        }
        res.myem<-myem.mix(dataset = data.na.intern,
                           K = nb.clust,
                           silent = TRUE,
                           nbquali=length(quali.index))
        
        newtheta<-try(da.mix(res.myem$s,res.myem$theta,steps=ceiling(Lstart),prior=res.myem$theta$pi), silent = TRUE)
        if(inherits(newtheta,"try-error")){newtheta<-try(da.mix(res.myem$s,res.myem$theta,steps=ceiling(Lstart)),silent = TRUE)}
        for(tab in seq.int(m)){
          if(verbose){cat(tab,"...",sep="")}
          res.try<-try(da.mix(res.myem$s,newtheta,steps=L,prior=res.myem$theta$pi),silent = TRUE)
          if(inherits(res.try, "try-error")) {
            newtheta <- try(da.mix(res.myem$s, newtheta, steps = L), silent = TRUE)
            cond1 <- cond2 <- inherits(newtheta, "try-error")
            if (!cond1) {
              cond2 <- any(unlist(lapply(newtheta, is.nan)))
            }
            if (cond1 | cond2) {
              warning(
                paste0("Imputation using JM-GL fails to impute ", m," datasets. Try to change the seed argument or use another imputation method (e.g. FCS-homo)"))
              (break)()
            }
          }else if(any(unlist(lapply(res.try, is.nan)))){
            warning(
              paste0("Imputation using JM-GL fails to impute ", m," datasets. Try to change the seed argument or use another imputation method (e.g. FCS-homo)"))
            (break)()
          }else{
            newtheta <- res.try
          }
          res.imp[[tab]]<-as.data.frame(imp.mix(res.myem$s,newtheta))
          #on stocke la partition
          temp.part[tab,] <- res.imp[[tab]]$class
          
          res.imp[[tab]]$class<-NULL
          if(data.type=="mixed"){
            res.imp[[tab]][,names(quali.index)]<-lapply(res.imp[[tab]][,names(quali.index),drop=FALSE],as.factor)
            for(qualitmp in names(quali.index)){
              levels(res.imp[[tab]]$qualitmp)<-levels(data.na$qualitmp)
            }
            res.imp[[tab]]<-res.imp[[tab]][,colnames(data.na)]
          }
        }
        maxm <- sum(!sapply(res.imp,is.null))
        res.conv<-sapply(seq.int(maxm),FUN = function(ii,res.imp,part,quanti.index){
          apply(res.imp[[ii]][,quanti.index,drop=FALSE],2,
                FUN=function(xx,ref){
                  res.by<-do.call(rbind,by(xx,INDICES=ref,
                                           FUN=function(yy){
                                             meangp<-mean(yy)
                                             fregp<-length(yy)
                                             res.out<-c(mean=meangp,freq=fregp)
                                             return(res.out)
                                           }))
                  inter<-sum(((res.by[,"mean"]-mean(xx))^2)*(res.by[,"freq"]))/sum(res.by[,"freq"])
                  return(inter)
                },ref=part[ii,])
        },res.imp=res.imp,
        part=temp.part,
        quanti.index= quanti.index) 
    }else{
      data.na.intern<-sapply(data.na,function(xx){as.numeric(xx)})
      res.myem<-myem.cat(dataset = data.na.intern,
                         K = nb.clust,
                         silent = TRUE)
      newtheta<-try(da.cat(res.myem$s,res.myem$theta,steps=Lstart,showits = (verbose==2)),silent=TRUE)
      for(tab in seq.int(m)){
        if(verbose){cat(tab,"...",sep="")}
        newtheta<-try(da.cat(res.myem$s,newtheta,steps=ceiling(L),showits = (verbose==2)),silent=TRUE)

        res.imp[[tab]]<-as.data.frame(imp.cat(res.myem$s,newtheta))
        res.imp[[tab]]$class<-NULL
        
          res.imp[[tab]]<-lapply(res.imp[[tab]],as.factor)
          for(qualitmp in colnames(data.na)){
            levels(res.imp[[tab]]$qualitmp)<-levels(data.na$qualitmp)
          }
          res.imp[[tab]]<-as.data.frame(res.imp[[tab]])
        }
      
    }
    if(verbose){cat("done!\n")}
  }
  if("JM-DP"== method){
    if(data.type=="continuous"){
      if(L==1){warning("The current version of JM-DP requires L>1")}
    data_obj <- readData(Y_in = data.na,99)
    model_obj <- createModel(data_obj, K_mix_comp = nb.clust)
    result_obj <- multipleImp(model_obj = model_obj, data_obj=data_obj,n_burnin = Lstart, m_Imp =m, interval_btw_Imp = L,show_iter = verbose)
    res.conv<-sapply(seq.int(m),FUN = function(ii,res.imp,part){
      apply(res.imp[ii,,],2,
            FUN=function(xx,ref){
              res.by<-do.call(rbind,by(xx,INDICES=ref,
                                       FUN=function(yy){
                                         meangp<-mean(yy)
                                         fregp<-length(yy)
                                         res.out<-c(mean=meangp,freq=fregp)
                                         return(res.out)
                                       }))
              inter<-sum(((res.by[,"mean"]-mean(xx))^2)*(res.by[,"freq"]))/sum(res.by[,"freq"])
              return(inter)
            },ref=part[ii,])
    },res.imp=result_obj$multiple_Imp,part=result_obj$multiple_Z_vec)
    
    rownames(res.conv)<-colnames(data.na)
    res.imp<-list()
    for(tab in seq.int(m)){
      res.imp[[tab]]<-as.data.frame(result_obj$multiple_Imp[tab,,])
      colnames(res.imp[[tab]])<-colnames(data.na)
    }
    }else if(data.type=="categorical"){
      res.imp<-DPMPM_nozeros_imp(X = data.na,nrun = Lstart+L*m,burn = Lstart,
                        thin = 1,K = nb.clust,aalpha = 0.25,
                        balpha = 0.25,m = m,seed = seed,silent = !verbose)$impdata
    }else if(data.type=="mixed"){
      stop("Not available. The MixedDataImpute package was removed from the CRAN repository. Please use another imputation method.")
        }
    if(verbose){cat("done!\n")}
  }
  res.out<-list(res.imp=res.imp,res.conv=res.conv,call=Call)
  return(res.out)
}