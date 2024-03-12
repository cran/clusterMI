#' @importFrom knockoff create.fixed stat.glmnet_coefdiff knockoff.filter
#' @importFrom withr with_seed
#' @importFrom mix rngseed prelim.mix imp.mix em.mix
#' @importFrom glmnet cv.glmnet

  varselbest.intern<-function(nnodes,
                            X,#une matrice
                            Y,#un vecteur
                            B=3000,
                            path.outfile=NULL,#chemin pour affichage
                            methods="knockoff",#c("glmnet","stepwise","knockoff")
                            sizeblock=ceiling(nrow(X)/10),#nb de variable choisies par simu
                            printflag=TRUE,#un peu d'affichage
                            r=NULL,#seuil une fois les B simu effectu?es
                            seed=1234567,
                            knockoff.arg=list(knockoffs=create.fixed,
                                              statistic=stat.glmnet_coefdiff,
                                              fdr=0.05,
                                              offset=0),
                            glmnet.arg=list(type.measure='mse',nfolds.cv=5,
                                            family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),
                                            offset=NULL, alpha = 1, nlambda = 100,
                                            lambda.min.ratio = ifelse(nrow(X)<(ncol(X)+1),0.01,0.0001), lambda=NULL,
                                            standardize = TRUE, intercept=TRUE, thresh = 1e-07,  dfmax = (ncol(X)+1) + 1,
                                            pmax = min((ncol(X)+1) + 1 * 2+20, (ncol(X)+1)), penalty.factor = rep(1, (ncol(X)+1)),
                                            lower.limits=-Inf, upper.limits=Inf, maxit=100000,
                                            type.gaussian=ifelse((ncol(X)+1)<500,"covariance","naive"),
                                            type.logistic=c("Newton","modified.Newton"),
                                            standardize.response=FALSE, type.multinomial=c("ungrouped","grouped")),
                            stepwise.arg=list(scale = 0,
                                              direction = c("both"),
                                              trace = 1),
                            nb.clust,modelNames
){
  # if(!setequal(order(sizeblock),seq(length(sizeblock)))){stop("les valeurs de sizeblock doivent etre decroissantes")}
 
  

  
  # tmp<-list(...)
  # if("kpmm"%in%names(tmp)){kpmm<-tmp[["kpmm"]]}else{kpmm<-5}
  # 
  #generation varblock (ss-ensemble de variables sur lequel est effectu? la selection
  Call<-list("nnodes"=nnodes,
             "X"=X,#une matrice
             "Y"=Y,#un vecteur
             "B"=B,
             "path.outfile"=path.outfile,
             "methods"=methods,
             "sizeblock"=sizeblock,
             "printflag"=printflag,
             "r"=r,
             "seed"=seed,
             "knockoff.arg"=knockoff.arg,
             "glmnet.arg"=glmnet.arg,
             "stepwise.arg"=stepwise.arg,
             "nb.clust"=nb.clust,
             "modelNames"=modelNames)
  if(!is.null(seed)){set.seed(seed)}
  if(is.data.frame(X)){X.intern<-as.matrix(X)}else{X.intern<-X}
  # print(dim(X.intern))
  # print(ceiling(ncol(X.intern)/sizeblock))
  # print(sizeblock)
  B.tmp<-0
  listvarblock<-list()
  while(B.tmp<B){
    
    #on repartit les variables dans une matrice, chaque ligne d?finira le sous ensemble de variables pour une simu
    matvarblock<-matrix(NA,sizeblock,ceiling(ncol(X.intern)/sizeblock))
    # print(length(matvarblock))
    matvarblock[seq(ncol(X.intern))]<-sample(seq(ncol(X.intern)),
                                             size = min(ncol(X.intern),length(matvarblock)),
                                             replace=FALSE)
    matvarblock<-t(matvarblock)
    # print( matvarblock)
    if(any(is.na(matvarblock))){
      # warning("the last block has been omited because of non-convenient size")
      matvarblock<-na.omit(matvarblock)
    }
    
    B.tmp<-B.tmp+nrow(matvarblock)
    listvarblock<-c(listvarblock,by(matvarblock,
                                    INDICES = seq(nrow(matvarblock)),
                                    FUN=function(x){x}))
  }
  #chaque ?l?ments de lise varblock est un sous ensemble de variables avec un minimum de chevauchements
  # print(names(listvarblock))
  
  cl <- parallel::makeCluster(nnodes, type = "PSOCK")
  parallel::clusterExport(cl, list("algo.intern","imputation.intern","lm","step",
  "X.intern","Y","methods","listvarblock", "path.outfile","printflag", "knockoff.arg", "glmnet.arg","stepwise.arg","seed","nb.clust",
  "create.fixed", "stat.glmnet_coefdiff", "knockoff.filter","with_seed", "rngseed", "prelim.mix", "imp.mix", "em.mix", "cv.glmnet","modelNames"
  ), envir = environment())
 
  # parallel::clusterExport(cl, list("algo.intern","imputation.intern","lm","step",
                                   # "X.intern","Y","methods","listvarblock", "path.outfile","printflag", "knockoff.arg", "glmnet.arg","stepwise.arg","seed","nb.clust"), envir = environment())
  #parallel::clusterEvalQ(cl, library("mix"));parallel::clusterEvalQ(cl, library("glmnet"));parallel::clusterEvalQ(cl, library("mice"));parallel::clusterEvalQ(cl, library("knockoff"));parallel::clusterEvalQ(cl, library("withr"))
  if (!is.null(path.outfile)) {
    parallel::clusterEvalQ(cl, sink(paste0(path.outfile, "/output", 
                                           Sys.getpid(), ".txt")))
  }
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, seed)
    parallel::clusterEvalQ(cl,rngseed(seed))  #set random number generator seed
  }
  
  res.par <- parallel::parSapply(cl, as.list(1:B), FUN = function(ii,
                                                                  XX,
                                                                  Y,
                                                                  methods,
                                                                  listvarblock,
                                                                  printflag,
                                                                  knockoff.arg,
                                                                  glmnet.arg,
                                                                  stepwise.arg,
                                                                  nb.clust,modelNames){
    
    if(printflag){cat(ii,"...")}
    
    res <- algo.intern(Xtrou=XX,Y,methods,varblock=as.matrix(listvarblock[[ii]]),
                    knockoff.arg=knockoff.arg,
                    glmnet.arg=glmnet.arg,
                    stepwise.arg=stepwise.arg,
                    nb.clust=nb.clust, modelNames=modelNames)
  },XX=X.intern, Y=Y, methods=methods, listvarblock=listvarblock, printflag=printflag, knockoff.arg=knockoff.arg, glmnet.arg=glmnet.arg, stepwise.arg=stepwise.arg, nb.clust=nb.clust, modelNames=modelNames, simplify = FALSE)
  
  parallel::stopCluster(cl)
  
  gardeprop<-lapply(lapply(res.par,"[[","res.select"),FUN="[[","garde")
  tmp<-do.call(rbind,gardeprop)
  # tmp[tmp<0]<-0
  garde<-colSums(tmp)
  # garde<-colSums(do.call(rbind,gardeprop))
  
  
  selectprop<-lapply(lapply(res.par,"[[","res.select"),FUN="[[","SousSelect")
  SousSelect<-colSums(do.call(rbind,selectprop))
  
  failureprop<-lapply(lapply(res.par,"[[","res.select"),FUN="[[","failure")
  failure<-colSums(do.call(rbind,failureprop))
  
  
  proportion<-SousSelect/garde
  names(proportion)<-names(garde)
  
  if(is.null(r)){selection<-NULL}else{selection<-colnames(X)[proportion>=r]}
  
  res<-list(garde=garde,
            effectif=SousSelect,
            proportion=proportion,
            selection=selection,
            failure=failure,
            listvarblock=listvarblock)
  if(printflag){cat("done!\n")}
  return(list(res=res,res.detail=res.par,call=Call))
}