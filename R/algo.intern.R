#' @importFrom withr with_seed
#' @importFrom mix rngseed
#' @importFrom knockoff knockoff.filter
#' @importFrom glmnet cv.glmnet
#' @importFrom mix prelim.mix imp.mix em.mix
#' @importFrom mclust hc Mclust

algo.intern<-function(Xtrou,Y,
                      methods,
                      varblock,
                      knockoff.arg,
                      glmnet.arg,
                      stepwise.arg,
                      nb.clust,
                      modelNames){
  #prend en entree un sous ensemble de variables et effectue dessus la selectio selon une m?thode choisie
  Nvar<-ncol(Xtrou)
  SousSelect <- garde <- failure<- rep(0,length=Nvar)
  names(SousSelect) <- names(garde) <- names(failure)<-colnames(Xtrou)
  if(sum(is.na(Xtrou))==0){method.na.tmp<-"cc"}else{method.na.tmp<-"mixture"}#pour g?rer le cas complet
  sousX = Xtrou[,as.vector(varblock)]
  if(method.na.tmp=="cc"){
    qui <- apply(!is.na(sousX),1,all)
    sousX <- sousX[qui,,drop=FALSE] ## et seulement les observations compl?tes
    res.init.hc <- hc(sousX,
                      modelName = "VVV",
                      use = "STD")
    clustering <- Mclust(sousX,
                         G=nb.clust,
                         initialization = list(hcPairs = res.init.hc),
                         modelNames=modelNames)$classification
    if(methods=="knockoff"){if (nrow(sousX)<=2*ncol(sousX)) return(NULL)}
    if(methods!="knockoff"){if (nrow(sousX)<=ncol(sousX)) return(NULL)}
  }else{
    qui <- seq(nrow(sousX))
    res.imp.intern<-imputation.intern(cbind(sousX,Y),nb.clust = nb.clust)
    sousX <-res.imp.intern$ximp[,seq(ncol(sousX))]
    clustering<-res.imp.intern$cluster
  }
  
  garde[varblock] <- garde[varblock]+1
  #centrage par groupe
  meanpergp<-unlist(by(Y[qui],INDICES=clustering,FUN = mean,simplify = FALSE,na.rm=TRUE))
  meanpergp[is.na(meanpergp)]<-mean(Y[qui],na.rm=TRUE)#si pas de valeur de Y dans un groupe, on prend la moyenne globale
  meanperind<-meanpergp[clustering]  
  sousY<-Y[qui]-meanperind

  ##########
  #knockoff
  ##########
  if("knockoff"==methods[1]){
    tmp<-try(with_seed(0,knockoff.filter(cbind(cluster=clustering,sousX),
                                                          sousY,
                                                          fdr = knockoff.arg$fdr,
                                                          offset=knockoff.arg$offset,
                                                          knockoffs=knockoff.arg$knockoffs,
                                                          statistic=knockoff.arg$statistic)),silent=FALSE)

    
    # tmp<-try(with_seed(0,knockoff.filter(sousX,sousY,
    # fdr = 0.1,
    # offset=0,
    # knockoffs=create.fixed.equi,
    # statistic=stat.lasso_signed_max)),silent=FALSE)
    
    if(!inherits(tmp,"try-error")){
      tmpout = tmp$selected
      SousSelect[which(colnames(Xtrou)%in%names(tmpout))] <- SousSelect[which(colnames(Xtrou)%in%names(tmpout))]+1
    }else{
      garde[varblock] <- garde[varblock]-1
      failure[varblock] <- failure[varblock]+1
    }
    res.select<-list(SousSelect=SousSelect, garde=garde, failure=failure)
  }
  
  
  
  ###########
  #SousLasso
  ##########
  
  if("glmnet"==methods[1]){
    glmnet1 = try(cv.glmnet(cbind(clustering,sousX),
                                    sousY,
                                    type.measure=glmnet.arg$type.measure.cv,
                                    nfolds=glmnet.arg$nfolds.cv,
                                    family=glmnet.arg$family,
                                    offset=glmnet.arg$offset, alpha = glmnet.arg$alpha, nlambda = glmnet.arg$nlambda,
                                    lambda.min.ratio = glmnet.arg$lambda.min.ratio, lambda=glmnet.arg$lambda,
                                    standardize = glmnet.arg$standardize, intercept=glmnet.arg$intercept, thresh = glmnet.arg$thresh,  dfmax = glmnet.arg$dfmax,
                                    pmax = glmnet.arg$pmax, penalty.factor = glmnet.arg$penalty.factor,
                                    lower.limits=glmnet.arg$lower.limits, upper.limits=glmnet.arg$upper.limits, maxit=glmnet.arg$maxit,
                                    type.gaussian=glmnet.arg$type.gaussian,
                                    type.logistic=glmnet.arg$type.logistic,
                                    standardize.response=glmnet.arg$standardize.response, type.multinomial=glmnet.arg$type.multinomial))#modif ici alpha (qui etait specifi?) est un parametre pour l'eslasticnet + probleme dans la fonction cv.elnet observ? qd k=2
    if(!inherits(glmnet1,"try-error")){
      cc <- as.matrix(coef(glmnet1,s='lambda.min',exact=TRUE))
      tmp <- row.names(cc)[which(cc!=0)]
      tmp <- which(colnames(Xtrou)%in%tmp)
      SousSelect[tmp] <- SousSelect[tmp]+1
    }else{
      garde[varblock] <- garde[varblock]-1
      failure[varblock] <- failure[varblock]+1
    }
    res.select<-list(SousSelect=SousSelect, garde=garde, failure=failure)
  }
  
  #########
  #stepwise
  #########
  if("stepwise"==methods[1]){
    sousX<-cbind(clustering,as.matrix(sousX))
    tmp <- try(step(lm(sousY~sousX),
                    trace= stepwise.arg$trace,
                    direction= stepwise.arg$direction,
                    scale=stepwise.arg$scale)) ##Y et les variaables selectionnees par tmp
    if(!inherits(tmp,"try-error")){  
      tmp <- names(tmp$coefficients)[-c(1,2)] ##-1 pour intercept
      tmp <- as.numeric(substr(tmp,7,10))
      SousSelect[tmp] <- SousSelect[tmp]+1
    }else{
      garde[varblock] <- garde[varblock]-1#on annule la simu
      failure[varblock] <- failure[varblock]+1
    }
    res.select<-list(SousSelect=SousSelect, garde=garde, failure=failure)
  }
  
  res<-list(res.select=res.select,qui=qui)
  return(res)
}