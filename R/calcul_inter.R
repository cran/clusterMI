#' @keywords internal
calcul_inter<- function(matrixpart){
  # matrixpart : matrice des partitions positionnées en colonnes (n x M)
  distance_from_idxs <- function(idxs,matrixpart) {
    # idxs vecteur de taille 2 correspondant à deux indice de colonnes cd matrixpart
    i1 <- idxs[1]
    i2 <- idxs[2]
    #calcul randindex
    v1 <- as.vector(matrixpart[,i1 ])#part 1
    v2 <- as.vector(matrixpart[,i2 ])#part 2
    tab <- table(v1,v2)
    cond1 <- all(dim(tab) == c(1, 1))
    if (cond1) {
      ari<- 1
    }else if(!cond1){
      a <- sum(choose(tab, 2))
      b <- sum(choose(rowSums(tab), 2)) - a
      c <- sum(choose(colSums(tab), 2)) - a
      d <- choose(sum(tab), 2) - a - b - c
      ari <- (a + d)/(a + b + c + d)
    }
    #calcul de la dissimilarité
    res.out <- (1 - ari)
    return(res.out)
  }
  size <- ncol(matrixpart)
  d <- apply(utils::combn(size, 2), 2, distance_from_idxs,matrixpart=matrixpart)
  varinter<-mean(d)
  return(varinter)
}
