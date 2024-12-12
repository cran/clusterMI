#' internal function
#'
#' @description prend en entrée un jeu de données, ajoute une variable classe, et renvoie les parametre du general location model. Une des particularité est l'initialisation du EM par un modèle de mélange. En effet em.mix a sa propre façon d'initialiser qui parfois renvoie des paramètres identiques pour plusieurs groupes
#' @importFrom mix prelim.mix em.mix
#' @importFrom mclust Mclust
#' @keywords internal

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
  
  cond1 <- (ncol(dataset) > nbquali) & (nbcompleteobs > length(quanti.index))
  cond2 <- (ncol(dataset) > nbquali) & (nbcompleteobs <= length(quanti.index))
  if (cond1 | cond2) {
    nbmod <- 1
    if (nbquali > 0) {
      nbmod <- prod(apply(dataset[, seq.int(nbquali), 
                                  drop = FALSE], 2, FUN = function(xx) {
                                    length(table(xx))
                                  }))
    }
  }
  if (cond1) {
    res.mclust <- Mclust(na.omit(dataset[, quanti.index, 
                                         drop = FALSE]),
                         G = K * nbmod,
                         modelNames = "EEI", 
                         verbose = FALSE)
    #transformation de la sortie de mclust pour s'insérer dans em.mix (réciproque de getparam.mix)
    start <- vector("list", length = 3)
    names(start) <- c("sigma","mu","pi")
    start$pi <- res.mclust$parameters$pro
    start$mu <- (res.mclust$parameters$mean - matrix(s$xbar, 
                                                     ncol = K * nbmod, nrow = length(quanti.index)))/matrix(s$sdv, 
                                                                                                            ncol = K * nbmod, nrow = length(quanti.index))
    mattmp <- res.mclust$parameters$variance$Sigma
    prodsd<-(matrix(s$sdv, s$q, s$q))*t((matrix(s$sdv, s$q, s$q)))
    mattmp <- mattmp/prodsd
    start$sigma <- mattmp[lower.tri(mattmp, diag = TRUE)]
    thetahat <- try(em.mix(s, showits = !silent, eps = eps, 
                           start = start),silent=TRUE)
    
    if (inherits(thetahat, "try-error")) {
      thetahat <- em.mix(s, showits = !silent, eps = eps)
    }
  }
  else if (cond2) {
    continue <- TRUE
    timestart<-Sys.time()
    timemax <- 15
    while(continue){
    #imputation de dataset par tirage aléatoire
    don.imp <- dataset[, quanti.index, drop = FALSE]
    valuesimp <- as.vector(unlist(sapply(as.data.frame(don.imp), 
                                         FUN = function(xx) {
                                           sample(na.omit(xx), size = sum(is.na(xx)), replace = TRUE)
                                         })))
    don.imp[is.na(don.imp)] <- valuesimp
    res.mclust <- Mclust(don.imp, G = K * nbmod, modelNames = "EEI", 
                         verbose = FALSE)
    #transformation de la sortie de mclust pour s'insérer dans em.mix (réciproque de getparam.mix)
    start <- vector("list", length = 3)
    names(start) <- c("sigma","mu","pi")
    start$pi <- res.mclust$parameters$pro
    start$mu <- (res.mclust$parameters$mean - matrix(s$xbar, 
                                                     ncol = K * nbmod, nrow = length(quanti.index)))/matrix(s$sdv, 
                                                                                                            ncol = K * nbmod, nrow = length(quanti.index))
    mattmp <- res.mclust$parameters$variance$Sigma
    prodsd<-(matrix(s$sdv, s$q, s$q))*t((matrix(s$sdv, s$q, s$q)))
    mattmp <- mattmp/prodsd
    start$sigma <- mattmp[lower.tri(mattmp, diag = TRUE)]
    thetahat <- try(em.mix(s, showits = !silent, eps = eps, 
                           start = start),silent=TRUE)
    continue <- inherits(thetahat, "try-error")&(difftime(Sys.time(),timestart)<=timemax)
    }
  }  else {
    thetahat <- em.mix(s, showits = !silent, eps = eps)
  }
  res.out <- list(s = s, theta = thetahat, ind.boot = ind.boot)
  return(res.out)
}
