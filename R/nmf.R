#' Consensus clustering using non-negative matrix factorization
#' 
#' @param listpart a list of partitions
#' @param nb.clust an integer specifying the number of clusters
#' @param threshold a real specifying when the algorithm is stoped. Default value is 10^(-5)
#' @param printflag a boolean. If TRUE, nmf will print messages on console. Default value is TRUE
#' @param nstart how many random sets should be chosen for kmeans initalization. Default value is 100
#' @param iter.max the maximum number of iterations allowed for kmeans. Default value is 50
#' @seealso \code{\link[stats]{kmeans}}
#' @keywords internal

nmf<-function(listpart,nb.clust,threshold=10^(-5),printflag=TRUE,nstart=100,iter.max=50){
  # print(str(listpart))
  listpart.intern<-listpart[!sapply(listpart,is.null)]
if(length(listpart.intern)==1){
  res<-list(Htilde=NULL,S=NULL,crit=NULL,Mtilde=NULL,cluster=listpart.intern[[1]])
  return(res)
}else if(length(listpart.intern)==0){ res<-list(Htilde=NULL,S=NULL,crit=NULL,Mtilde=NULL,cluster=NULL);return(res)}
  
  Mtilde<-matrix(0,length(listpart.intern[[1]]),length(listpart.intern[[1]]))#matrice d'adjacence rempli de 0
  for(ii in seq(length(listpart.intern))){
    Mtilde<-Mtilde+tcrossprod(FactoMineR::tab.disjonctif(listpart.intern[[ii]]))#somme de toutes les matrices d'adjacences
  }
  Mtilde<-Mtilde/length(listpart.intern)#moyenne des matrices d'adjacences

  #initialisation
  res.kmeans<-try(kmeans(Mtilde,centers = nb.clust,nstart = nstart,iter.max=iter.max))
if("try-error"%in% class(res.kmeans)){
  res.kmeans<-sample(seq(nb.clust),size = ncol(Mtilde),replace=TRUE)
}else{
    res.kmeans<-res.kmeans$cluster
}

  H<-FactoMineR::tab.disjonctif(res.kmeans)
  S<-crossprod(H)
  Htilde<-H%*%diag(diag(1/sqrt(S)),nb.clust,nb.clust)
  continue<-TRUE
  critsave<-sqrt(sum(diag(crossprod(Mtilde-Htilde%*%tcrossprod(S,Htilde)))))
  comp<-1
  while(continue){
    if(printflag){cat(comp,"...")}
  #Htilde update
  MHS <- Mtilde%*%Htilde%*%S
  multHtilde<-sqrt(MHS/(tcrossprod(Htilde)%*%MHS))
  multHtilde[is.nan(multHtilde)]<-0
  Htilde<-Htilde*multHtilde
  #S update
  cpHtilde <- crossprod(Htilde)
  multS<-sqrt((crossprod(Htilde,Mtilde)%*%Htilde)/(cpHtilde%*%S%*%cpHtilde))
  multS[is.nan(multS)]<-0
  multS[is.infinite(multS)]<-0
  S<-S*multS
  
  critsave<-c(critsave,sqrt(sum(diag(crossprod(Mtilde-Htilde%*%tcrossprod(S,Htilde))))))
  comp<-comp+1
  diffcrit<-critsave[comp-1]-critsave[comp]
  continue<-(diffcrit>=threshold)
  }
  if(printflag){cat("done \n")}
  
 res<-list(Htilde=Htilde,S=S,crit=critsave,Mtilde=Mtilde,cluster=apply(Htilde,1,which.max))
 return(res)
}