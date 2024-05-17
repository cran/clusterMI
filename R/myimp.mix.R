#' @importFrom mix prelim.mix em.mix imp.mix rngseed
#' @import mix
#' @importFrom mclust Mclust

myimp.mix<-function(dataset,K,seed=1234567,boot=FALSE,silent=TRUE,nbquali=0, eps=0.0001){
  # if(silent){sink("tmp.txt")}
  nbcompleteobs <- sum(apply(!is.na(dataset), 1, all))
  quanti.index <- seq.int(from = (nbquali + 1), to = ncol(dataset))
  
  if(is.data.frame(dataset)){
    don<-as.matrix(dataset)
  }else{don<-dataset}
  
  donclass<-cbind(class=c(K,rep(NA,nrow(don)-1)),don)
  s <- prelim.mix(donclass,1+nbquali)
  if(!boot){
    if ((ncol(dataset) > nbquali)& (nbcompleteobs > length(quanti.index))) {
      #si presence de variables quanti et suffisament d'individus complets
      
      #calcul du nombre de modalités
      nbmod <- 1
      if (nbquali > 0) {
        nbmod <- prod(apply(dataset[, seq.int(nbquali), 
                                    drop = FALSE], 2, FUN = function(xx) {
                                      length(table(xx))
                                    }))
      }
      #initialisation par modèle de mélange (permet d'éviter parametre identique entre différents groupes dans em.mix)
      res.mclust <- Mclust(na.omit(dataset[, quanti.index, 
                                           drop = FALSE]), G = K * nbmod, modelNames = "EEI", 
                           verbose = FALSE)
      #transformation de la sortie de Mclust selon la réciproque de mix:getparameter
      start <- list()
      start$pi <- res.mclust$parameters$pro
      start$mu <- (res.mclust$parameters$mean - matrix(s$xbar, 
                                                       ncol = K * nbmod, nrow = length(quanti.index)))/matrix(s$sdv, 
                                                                                                              ncol = K * nbmod, nrow = length(quanti.index))
      mattmp <- res.mclust$parameters$variance$Sigma/t(matrix(s$sdv, 
                                                              s$q, s$q))/matrix(s$sdv, s$q, s$q)
      start$sigma <- mattmp[upper.tri(mattmp, diag = TRUE)]
      #EM
      thetahat <- try(em.mix(s, showits = !silent, eps = eps, 
                             start = start),silent=TRUE)
      if(inherits(thetahat,"try-error")){
        thetahat <- em.mix(s, showits = !silent, eps = eps)
      }
    }else {
      thetahat <- em.mix(s, showits = !silent, eps = eps)
    }
    ##############################
    # #initialisation pour em
    #   res.mclust <- Mclust(don, G = K, modelNames = "EEI", 
    #                        verbose = FALSE)
    #   start <- list()
    #   start$pi <- res.mclust$parameters$pro
    #   start$mu <- (res.mclust$parameters$mean - matrix(s$xbar, 
    #                                                    ncol = K, nrow = ncol(don)))/matrix(s$sdv, ncol = K, 
    #                                                                                        nrow = ncol(don))
    #   mattmp <- res.mclust$parameters$variance$Sigma/t(matrix(s$sdv, 
    #                                                           s$q, s$q))/matrix(s$sdv, s$q, s$q)
    #   start$sigma <- mattmp[upper.tri(mattmp, diag = TRUE)]
    #   #em
    #   thetahat <- em.mix(s, showits = !silent, eps = eps, 
    #                      start = start)
      ##############################
  }else{
    ind.boot<-sample(seq(nrow(don)),
                     size=nrow(don),
                     replace=TRUE)
    donclassboot<-cbind(class=c(K,rep(NA,nrow(don)-1)),
                        don[ind.boot,])
    
    sboot <- prelim.mix(donclassboot,1+nbquali)
    # s$z<-donclass
    thetahat <- em.mix(sboot,showits = !silent, eps=eps)
  }
  
  
  if(!is.null(seed)){mix::rngseed(seed)}     # set random number generator seed
  ximp <- imp.mix(s, thetahat, donclass)
  # if(silent){sink();file.remove("tmp.txt")}
  if(!boot){
    return(ximp)
  }else{
    res.out<-list(ximp=ximp,theta=thetahat,s=s)
    return(res.out)
  }
}