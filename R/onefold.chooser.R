#' one fold cross-validation for specifying threshold r
#' @importFrom stats lm predict.lm
#' @keywords internal
#' @param data.train.k a train set
#' @param data.test.k a test set
#' @param jj name of the outcome used for variable selection
#' @param grid.r a grid for the tuning parameter r
#' @param nb.clust number of clusters
#' @param sizeblock number of sampled variables at each iteration
#' @param method.select variable selection method
#' @param B number of iterations
#' @param modelNames mixture model specification for imputation of subsets
#' @param K number of fold
#' @param path.outfile a path for message redirection
#' @param nbvarused a maximal number of selected variables (can be required with a large number of variables) 

onefold.chooser<-function(data.train.k,
                          data.test.k,
                          jj,#nom de la cible
                          grid.r,
                          nb.clust,
                          nnodes,
                          sizeblock,
                          method.select,
                          B,
                          modelNames,
                          K,
                          path.outfile,
                          nbvarused){
  rmin<- -Inf #pour affichage
  res.varsel.k<-varselbest(data.na=data.train.k,
                           nb.clust=nb.clust,
                           listvar = jj,
                           nnodes = nnodes,
                           sizeblock = sizeblock,
                           method.select = method.select,
                           B=B,
                           graph = FALSE,
                           printflag = FALSE,
                           r=0,
                           modelNames = modelNames,
                           path.outfile = path.outfile)
  #gestion gde dim

  #on ajuste la grille en fonction du nombre de variables
  #on fait coïncider les proportions et les valeurs de la grille en arrondissant
  nbdec.grid<-nchar(as.character((length(grid.r)-1)))
  roundprop<-round(res.varsel.k$proportion[jj,],digits = nbdec.grid)
  grid.r.intern.k<-grid.r[grid.r>=quantile(roundprop,
                                           1-(min(1,#on prend le quantile minimum de la grille
                                                  (((nrow(data.train.k)+nrow(data.test.k))*((K-1)/K))+1)/((ncol(data.train.k)-1)),#on prend le quantile minimum pour gérer la gde dim
                                                  nbvarused/(ncol(data.train.k)-1)#on prend le quantile pour satisfaire le nombre choisi
                                           )
                                           ))]
  
  #gestion NA (en fusionnant)
  
  data.train.k.small<-data.train.k[,which((res.varsel.k$proportion[jj,]>=grid.r.intern.k[1])&(!(colnames(res.varsel.k$proportion)==jj))),drop=FALSE]
  data.test.k.small<-data.test.k[,which((res.varsel.k$proportion[jj,]>=grid.r.intern.k[1])&(!(colnames(res.varsel.k$proportion)==jj))),drop=FALSE]
  X.intern<-rbind(data.train.k.small,data.test.k.small)
  
  if(any(is.na(X.intern))){
    #imputation de xtrain, Xtest et Y train simultanement / ajustement / prevision
    X.tmp <- cbind(X.intern,Y=c(data.train.k[,jj],rep(NA,nrow(data.test.k.small))))
    # save(X.tmp,file="C:/Users/vince/Desktop/tmp/Xtmp.Rdata")
    res.imputation.intern <- imputation.intern(X.tmp,nb.clust=nb.clust,showits=FALSE)
    if(is.null(res.imputation.intern)){
      warning("The EM algorithm does not converge. You should probably increase K.")
      res.imputation.intern <- imputation.intern(X.tmp,nb.clust=nb.clust,showits=FALSE,maxits=20)
    }
    clustering<- res.imputation.intern$cluster
    X.imp <- res.imputation.intern$ximp#imputation par JM-GL
    Xtrain.imp <- X.imp[seq(nrow(data.train.k.small)),seq(ncol(data.train.k.small)),drop=FALSE]
    Xtest.imp <- X.imp[-seq(nrow(data.train.k.small)),seq(ncol(data.train.k.small)),drop=FALSE]
  }else{
    res.init.hc <- hc(rbind(data.train.k.small,data.test.k.small),
                      modelName = "VVV",
                      use = "STD")
    clustering <- Mclust(rbind(data.train.k.small,data.test.k.small),
                         G=nb.clust,
                         initialization = list(hcPairs = res.init.hc),
                         modelNames=modelNames,verbose=FALSE)$classification
    Xtrain.imp<-data.train.k.small
    Xtest.imp<-data.test.k.small
  }
  clustering.train<-clustering[seq(nrow(data.train.k.small))]
  clustering.test<-clustering[-seq(nrow(data.train.k.small))]
  
  #on definit la liste des variables explicatives
  Xselectlist<-lapply(grid.r.intern.k,FUN=function(rr,propor){
    names(which(propor>=rr))
  },propor=res.varsel.k$proportion[jj,-which(colnames(res.varsel.k$proportion)==jj)])
  
  names(Xselectlist)<-grid.r.intern.k
  
  #on calcule l'erreur en prédiction sans les doublons
  sapply(unique(Xselectlist),length)
  errorunique<-colSums((sapply(unique(Xselectlist),
                               FUN = function(Xnames,Xtrain,Ytrain,Xtest,clustering.train,clustering.test){
                                 #ajustement
                                 if(length(Xnames)>0){
                                   res.lm<-try(lm(Y~.,data=data.frame("cluster"=as.factor(clustering.train),"Y"=Ytrain,Xtrain[,Xnames,drop=FALSE])))
                                   # print(res.lm)
                                   # res.lm<-try(lm(Ytrain~as.matrix(Xtrain[,Xnames])))
                                 }else{
                                   # res.lm<-try(lm(Ytrain~as.factor(clustering.train)))
                                   res.lm<-try(lm(Y~.,data=data.frame("cluster"=as.factor(clustering.train),"Y"=Ytrain)))
                                   
                                 }
                                 
                                 #prediction
                                 if("try-error"%in%class(res.lm)){
                                   return(rep(NA,nrow(Xtest)))
                                 }
                                
                                 res<-predict.lm(res.lm,
                                                 newdata = data.frame("cluster"=as.factor(clustering.test),
                                                                      Xtest[,Xnames,drop=FALSE]))
                                 
                                 return(res)
                               },
                               Xtrain=Xtrain.imp,
                               Ytrain=data.train.k[,jj],
                               Xtest=Xtest.imp,
                               clustering.test=clustering.test,
                               clustering.train=clustering.train)-matrix(data=rep(data.test.k[,jj],length(unique(Xselectlist))),
                                                       nrow=nrow(data.test.k),
                                                       ncol=length(unique(Xselectlist))))^2,na.rm=TRUE)
  # print(errorunique)
  # print(str(Xselectlist))
  # save(errorunique,file="C:/Users/vince/Desktop/tmp/errorunique.Rdata")
  # save(Xselectlist,file="C:/Users/vince/Desktop/tmp/Xselectlist.Rdata")
  #on duplique les erreurs pour chaque valeur de la grille
  error.k<-sapply(Xselectlist,
                  FUN = function(modelr,errorunique,nameerrorunique){
                    errorunique[which(sapply(nameerrorunique,FUN=setequal,y=modelr))]
                  },errorunique=errorunique,nameerrorunique=unique(Xselectlist),
                  simplify = FALSE)
  return(error.k)
 
}
