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
#'  \item{res.conv}{ for FCS methods, an array given the within inertia of each imputed variable at each iteration and for each imputed dataset}
#'  \item{call}{ the matching call}
#' 
#' @param data.na an incomplete dataframe
#' @param method a single string specifying the imputation method used among "FCS-homo","FCS-hetero","JM-DP","JM-GL". By default method = "JM-GL". See the details section
#' @param nb.clust number of clusters
#' @param m number of imputed datasets. By default, m = 20.
#' @param maxit number of iterations for FCS methods (only used for method = FCS-homo or method = FCS-hetero)
#' @param Lstart number of iterations for the burn-in period (only used if method ="JM-DP" or "JM-GL")
#' @param L number of skipped iterations to keep one imputed data set after the burn-in period (only used if method ="JM-DP" or "JM-GL")
#' @param method.mice a vector of strings (or a single string) giving the imputation method for each variable (only used for method = FCS-homo or method = FCS-hetero). Default value is "norm" for FCS-homo and "mice.impute.2l.jomo" for FCS-hetero
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
#' data(wine)
#' set.seed(123456)
#' wine.na <- wine
#' wine.na$cult <- NULL
#' wine.na <- as.matrix(wine.na)
#' wine.na[sample(seq(length(wine.na)), size = ceiling(length(wine.na)/3))] <- NA
#' nb.clust <- 3 # number of clusters
#' m <- 3 # number of imputed data sets
#' res.imp <- imputedata(data.na = wine.na, nb.clust = nb.clust, m = m)
#' lapply(res.imp$res.imp, summary)

imputedata<-function(data.na,method="JM-GL",
                    nb.clust=NULL,m=20,maxit=50,Lstart=100,L=20,method.mice=NULL,
                    predictmat=NULL,verbose=TRUE, seed=1234, bootstrap = FALSE){
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
  #check colnames (mice ne supporte pas n'importe quel nom)
  checkcolnames<-try(sapply(colnames(data.na),str2lang),silent=TRUE)
  if(("try-error"%in%class(checkcolnames))&method%in%c("FCS-homo","FCS-hetero")){
    stop("Colnames of data.na are not supported by mice. Please rename colnames with character string compatible with str2lang")
  }
  
  Call<-list(data.na=data.na,method=method,
  nb.clust=nb.clust,m=m,maxit=maxit,Lstart=Lstart,L=L,method.mice=method.mice,
  predictmat=predictmat,data.type=data.type,bootstrap=bootstrap)
  if(is.null(method.mice)){method.mice.homo<-NULL;method.mice.hetero<-"2l.jomo"}else{
  if("FCS-homo"%in%method){method.mice.homo<-method.mice}else if("FCS-hetero"%in%method){method.mice.hetero<-method.mice}}
  # if("FCS-homo"%in%method){bootstrap<-TRUE}else if("FCS-hetero"%in%method){bootstrap<-FALSE}
  # bootstrap<-FALSE
  res.conv<-NULL
  
  if(verbose){cat(method,"\n")}
  
  if(("FCS-homo"%in%method)){
    if(data.type=="continuous"){
      # data.na.intern<-data.na
      data.na.schaf<-data.na
      prior.dir<-try(myem.mix(dataset = data.na,
                          K = nb.clust,silent = TRUE)$theta$pi,silent=TRUE)#ML
      quali.index<-NULL
      res.conv<-array(NA,dim=c(m,maxit,ncol(data.na),2),dimnames=list(seq(m),seq(maxit),colnames(data.na),c("var_inter","var_intra")))
      
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
                      dimnames=list(seq(m),seq(maxit),colnames(data.na)[quanti.index],c("var_inter","var_intra")))
      
    }
    if("try-error"%in%class(prior.dir)){prior.dir<-0.5#;rngseed(seed)
    }
    res.imp<-list()
    for (tab in seq(m)) {
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
        if (("try-error" %in% class(data.init)) & (tab == 
                                                   1)) {
          mix::rngseed(seed)
        }
      }
      else if (data.type == "categorical") {
        data.init <- try(as.data.frame(myimp.cat(data.na.schaf, 
                                                 K = nb.clust, silent = TRUE, seed = seed.tmp)),silent = TRUE)
        if (("try-error" %in% class(data.init)) & (tab == 
                                                   1)) {
          rngseedcat(seed)
        }
      }
      if (!("try-error" %in% class(data.init))) {
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
        else {
          if (min(table(data.init$class)) <= min(rowSums(predictmat))) {
            flag <- TRUE
          }
        }
      }
      if (!("try-error" %in% class(data.init)) & (!flag)) {
        flag <- (length(table(data.init$class)) != nb.clust)
      }
      if (("try-error" %in% class(data.init)) | flag){
        cat("random initialization\n")
        data.init <- as.data.frame(data.na)
        data.init <- complete(mice(data = data.init, 
                                   maxit = 0))
        clustering.tmp <- sample(seq(nb.clust), size = nrow(data.na), 
                                 replace = TRUE)
        data.init.tmp <- cbind.data.frame(class = as.factor(clustering.tmp), 
                                          data.init)
      }else{
        clustering.tmp<-data.init$class
        data.init.tmp<-cbind.data.frame(class=as.factor(clustering.tmp),
                                        data.init[,-which(colnames(data.init)=="class")])
      }
      
      don.na.class<-cbind.data.frame(class=as.factor(clustering.tmp),data.na)
      if(verbose){cat(" nbiter = ")}
      for(nbiter in seq(maxit)){
        if(verbose){cat(nbiter,"...",sep="")}
        #besoin de mettre a jour le clustering a part, car completement manquant
        if(is.null(predictmat)){
          # on vérifie l'effectif des classes en fonction du nombre de variables
        if (min(table(don.na.class$class)) < ncol(don.na.class)) {
            warning("The number of individuals per cluster is too low to build the regression model. The imputation cycle is stopped.")
            break()
          }

          mice.tmp <- mice(data = don.na.class, m = 1, 
            method = method.mice.homo, maxit = 1, printFlag = FALSE, 
            data.init = data.init.tmp)
          if (!is.null((mice.tmp$loggedEvents))) {
            warning("mice returns a potential issue\n", mice.tmp$loggedEvents) 
          }

          data.init.tmp<-complete(mice.tmp)
        
        }else{
          # on vérifie l'effectif des classes en fonction du nombre de variables
          if (min(table(don.na.class$class)) <= min(rowSums(predictmat))) {
            warning("The number of individuals per cluster is too low to build the regression model. The imputation cycle is stopped.")
            break()
          }
          predictormatrix<-rbind(c(0,rep(0,nrow(predictmat))),cbind(rep(1,nrow(predictmat)),predictmat))
          colnames(predictormatrix)<-c("class",colnames(predictmat))
          rownames(predictormatrix)<-c("class",rownames(predictmat))
          data.init.tmp<-complete(mice(data = don.na.class,
                                       m = 1,
                                       method = method.mice.homo,
                                       maxit=1,
                                       printFlag=FALSE,
                                       data.init = data.init.tmp,predictorMatrix = predictormatrix))
        }
        
        varnotimputed<-names(which(apply(is.na(data.init.tmp),2,any)))
        if(length(varnotimputed)>0){warning(c("The following variables are not imputed: ",varnotimputed,". This can be due to colinearity between variables."))}
        if(data.type=="continuous"){
          data.init.tmp.intern<-data.init.tmp
          if(length(varnotimputed>0)){data.init.tmp.intern<-data.init.tmp[,-which(colnames(data.init.tmp)%in%varnotimputed)]}
        
        }else if(data.type=="categorical"){

          data.init.tmp.intern<-cbind.data.frame(class=data.init.tmp[,which(colnames(data.init.tmp)=="class")],
                                                 MCA(
                                                   data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))],
                                                   graph=FALSE,ncp=Inf)$ind$coord)
          
        }else if(data.type=="mixed"){
          data.init.tmp.intern<-cbind.data.frame(class=data.init.tmp[,which(colnames(data.init.tmp)=="class")],
                                                 FAMD(
                                                   data.init.tmp[,-which(colnames(data.init.tmp)%in%c("class",varnotimputed))],
                                                   graph=FALSE,ncp=Inf)$ind$coord)
          
        }
        if(!bootstrap){
        # tirage de W

        res.myem<-myem.mix(dataset = data.init.tmp.intern[,-which(colnames(data.init.tmp.intern)=="class")],
                           K = nb.clust,silent = TRUE)
        res.myimp<-myimp.mix(data.init.tmp.intern[,-which(colnames(data.init.tmp.intern)=="class")],
                             K = nb.clust,seed = NULL,silent = TRUE)
        W<-res.myimp[,"class"]
        # tirage de theta
        res.myem.new<-res.myem
        res.myem.new$s$w<-W
        newtheta<-da.mix(res.myem.new$s,res.myem.new$theta,steps=1,prior=prior.dir)
        res.myem.new$theta$pi<-newtheta$pi
        res.myem.new$s$w<-res.myem$s$w
        
        #imputation
        
        W<-imp.mix(res.myem.new$s,theta = res.myem.new$theta)[,"class"]
        
        }else{
          res.mclust<-mclustboot.intern(data.init.tmp.intern[,-which(colnames(data.init.tmp.intern)=="class")],
                                        G=nb.clust,modelNames = "EEE",verbose =FALSE)
          res.mclust$data<-as.matrix(data.init.tmp.intern[,-which(colnames(data.init.tmp.intern)=="class")])
          
          W<-drawW(res.mclust = res.mclust,method = "LDA")
          #gestion du cas ou une classe est vide
          if(length(unique(W))<nb.clust){
            tableW<-rep(0,nb.clust)
            names(tableW)<-as.character(seq(nb.clust))
            tableW[as.character(unique(W))]<-table(W)
          }else{
            tableW<-table(W)
          }
          
          
          theta<-as.vector(rdirichlet(1,alpha = prior.dir+tableW))
          res.mclust$parameters$pro<-theta
          W<-drawW(res.mclust = res.mclust,method = "LDA")
        }
        clustering.tmp<-W
        don.na.class<-cbind.data.frame(class=as.factor(as.character(W)),data.na)
        data.init.tmp$class<-as.factor(as.character(W))
        if(data.type%in%c("mixed","continuous")){
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
      }
      
      res.imp[[tab]]<-data.init.tmp
      res.imp[[tab]]$class<-NULL
    }
    if(verbose){cat("done!\n")}
  }
 
  if(("FCS-hetero"%in%method)){
    if(data.type!="continuous"){stop("FCS-hetero is not available for non-continuous data. Use FCS-homo (or a JM method).")}
    quali.index<-NULL
       prior.dir<-try(myem.mix(dataset = data.na,
                        K = nb.clust,silent = TRUE)$theta$pi,silent=TRUE)
    if("try-error"%in%class(prior.dir)){prior.dir<-0.5}
    res.conv<-array(NA,dim=c(m,maxit,ncol(data.na),2),dimnames=list(seq(m),seq(maxit),colnames(data.na),c("var_inter","var_intra")))
    res.imp<-list()
    # if(checkconv){par(mfrow=c(m,1))}
    for(tab in seq(m)){
      if(verbose){cat("\n m =",tab)}
      if(tab==1){data.init<-try(as.data.frame(myimp.mix(data.na,K = nb.clust,silent = TRUE, seed = seed)),silent=TRUE)
      }else{
        data.init<-try(as.data.frame(myimp.mix(data.na,K = nb.clust,seed=NULL,silent = TRUE)),silent=TRUE)
      }
      if("try-error"%in%class(data.init)){
        cat("random initialization\n")
        data.init<-as.data.frame(data.na)
        data.init<-cbind.data.frame(class=as.factor(sample(seq(nb.clust),size=nrow(data.na),replace=TRUE)),complete(mice(data = data.init,maxit=0)))
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
      for(nbiter in seq(maxit)){
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
          names(tableWstar)<-as.character(seq(nb.clust))
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
 
  if("JM-GL"%in%method){
    res.imp<-list()
    rngseed(seed)
    if(data.type%in%c("continuous","mixed")){
     
      data.na.intern<-data.na
        quali.index<-which(sapply(data.na,is.factor)|(sapply(data.na,is.ordered)))

        if(data.type=="mixed"){
          quanti.index<-which(!(sapply(data.na,is.factor)|(sapply(data.na,is.ordered))))
          data.na.intern[,quali.index]<-sapply(data.na.intern[,quali.index],as.numeric)
          data.na.intern<-data.na.intern[,c(quali.index,quanti.index)]
        }
        # print(data.na.intern)
        res.myem<-myem.mix(dataset = data.na.intern,
                           K = nb.clust,
                           silent = TRUE,
                           nbquali=length(quali.index))
        
        newtheta<-try(da.mix(res.myem$s,res.myem$theta,steps=ceiling(Lstart),prior=res.myem$theta$pi), silent = TRUE)
        if("try-error"%in%class(newtheta)){newtheta<-try(da.mix(res.myem$s,res.myem$theta,steps=ceiling(Lstart)),silent = TRUE)}
        for(tab in seq(m)){
          if(verbose){cat(tab,"...",sep="")}
          res.try<-try(da.mix(res.myem$s,newtheta,steps=L,prior=res.myem$theta$pi),silent = TRUE)
          if("try-error"%in%class(res.try)){
            newtheta<-try(da.mix(res.myem$s,newtheta,steps=L),silent=TRUE)
            if("try-error"%in%class(newtheta)){
              stop("Imputation using JM-GL fails. Try to change the seed argument or use another imputation method (i.e. FCS-homo)") 
              }
            }else{newtheta<-res.try}
          res.imp[[tab]]<-as.data.frame(imp.mix(res.myem$s,newtheta))
          res.imp[[tab]]$class<-NULL
          if(data.type=="mixed"){
            res.imp[[tab]][,names(quali.index)]<-lapply(res.imp[[tab]][,names(quali.index),drop=FALSE],as.factor)
            for(qualitmp in names(quali.index)){
              levels(res.imp[[tab]]$qualitmp)<-levels(data.na$qualitmp)
            }
            res.imp[[tab]]<-res.imp[[tab]][,colnames(data.na)]
          }
        }
    }else{
      data.na.intern<-sapply(data.na,function(xx){as.numeric(xx)})
      res.myem<-myem.cat(dataset = data.na.intern,
                         K = nb.clust,
                         silent = TRUE)
      newtheta<-try(da.cat(res.myem$s,res.myem$theta,steps=Lstart,showits = (verbose==2)),silent=TRUE)
      for(tab in seq(m)){
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
  if("JM-DP"%in%method){
    # require(DPImputeCont)
    if(data.type=="continuous"){
    data_obj <- readData(Y_in = data.na,99)
    model_obj <- createModel(data_obj, K_mix_comp = nb.clust)
    result_obj <- multipleImp(model_obj = model_obj, data_obj=data_obj,n_burnin = Lstart, m_Imp =m, interval_btw_Imp = L)
    res.imp<-list()
    for(tab in seq(m)){
      # if(verbose){cat(tab,"...",sep="")}
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