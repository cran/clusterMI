#' Compute Silhouette index
#' 
#' @param d a distance matrix
#' @param cl a vector of integer defining the partition
#' @return a real
#' @keywords internal

Silhouette.intern <- function(d,cl){
  #d matrice de distance (avec fonction "dist"
  #cl une partition
  d <- as.matrix(d)
  
  #gestion des groupes de taille 1
  Ai<-rep(NA,length(cl))
  classeone<-as.numeric(names(which(table(cl)==1)))
  indnotalone <- which(!(cl%in%classeone))
  Ai[cl%in%classeone] <- 1
  
  Ai[indnotalone] <- sapply(indnotalone,FUN = function(ii,d,cl){
    tmp<- d[ii,cl==cl[ii]]
    sum(tmp)/(length(tmp)-1)
    },d=d,cl=cl)
  
  Bi <- sapply(seq.int(length(cl)),FUN = function(ii,d,cl){
    j.grid<- unique(cl)
    j.grid<- j.grid[which(j.grid!=cl[ii])]
    min(sapply(j.grid,
               FUN = function(j,d,cl){
      mean(d[,cl==j])},cl=cl,d=d[ii,,drop=FALSE]))
  },d=d,cl=cl)
  
  res.out<-mean((Bi-Ai)/apply(cbind(Ai,Bi),1,max))
  return(res.out)
}
