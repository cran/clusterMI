#' @importFrom mix prelim.mix em.mix
#' @importFrom mclust Mclust

myem.mix<-function (dataset, K, boot = FALSE, silent = TRUE, nbquali = 0, 
          eps = 1e-04) 
{
  nbcompleteobs <- sum(apply(!is.na(dataset), 1, all))
  quanti.index <- seq.int(from = (nbquali + 1), to = ncol(dataset))
  ind.boot <- seq(nrow(dataset))
  if (is.data.frame(dataset)) {
    don <- as.matrix(dataset)
  }
  else {
    don <- dataset
  }
  don <- cbind(class = c(K, rep(NA, nrow(don) - 1)), don)
  if (!boot) {
    s <- prelim.mix(don, 1 + nbquali)
  }
  else {
    ind.boot <- sample(seq(nrow(don)), size = nrow(don), 
                       replace = TRUE)
    s <- prelim.mix(cbind(class = c(K, rep(NA, nrow(don) - 
                                             1)), don[ind.boot, ]), 1 + nbquali)
  }
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
  res.out <- list(s = s, theta = thetahat, ind.boot = ind.boot)
  return(res.out)
}
