#' @importFrom FactoMineR MCA FAMD
#' @importFrom fpc kmeansCBI nselectboot claraCBI noisemclustCBI hclustCBI
#' @importFrom stats dist

calculintra.intern<-function(res.imp,
                             method.clustering,
                             method.dist,Cboot,
                             nb.clust,scaling,
                             method.agnes,
                             nstart.kmeans,
                             modelNames,
                             m.cmeans){
  
  if(sum(sapply(res.imp,is.numeric))==ncol(res.imp)){
    #quanti
    res.imp.intern<-scale(res.imp,scale = scaling,center=FALSE)
  }else if(sum(sapply(res.imp,is.factor))==ncol(res.imp)){
    #quali
    res.imp.intern<-MCA(res.imp,ncp=Inf,graph=FALSE)$ind$coord
  }else{
    #mixte
    res.imp.intern<-FAMD(res.imp,ncp=Inf,graph=FALSE)$ind$coord
  }
  if(method.clustering=="kmeans"){
    clustermethod<-kmeansCBI
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=clustermethod,
                         classification="centroid",krange=nb.clust,scaling=FALSE,count = FALSE,runs=nstart.kmeans)$stabk[nb.clust]}
  if(method.clustering=="cmeans"){
    clustermethod<-cmeansCBI.intern
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=clustermethod,
                         classification="centroid",krange=nb.clust,scaling=FALSE,m.cmeans=m.cmeans)$stabk[nb.clust]}
  if(method.clustering%in%c("pam","clara")){
    clustermethod<-claraCBI
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=clustermethod,
                         classification="centroid",krange=nb.clust,usepam=(method.clustering=="pam"))$stabk[nb.clust]
  }
  if(method.clustering=="mixture"){
    if(is.null(modelNames)){
      methodclassif<- "qda"
    }else{
      methodclassif<-ifelse(test=modelNames%in%c("EEE","EEI","EII"),yes = "lda",no="qda")
    }
    
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=noisemclustCBI,
                         classification=methodclassif,krange=nb.clust,multipleboot=FALSE,modelNames=modelNames,verbose = FALSE)$stabk[nb.clust]
  }
  if(method.clustering=="hclust"){
    clustermethod<-hclustCBI
    res.out<-nselectboot(data=dist(res.imp.intern,method = method.dist),
                         B=Cboot,
                         distances = TRUE,
                         clustermethod=clustermethod,
                         classification="averagedist",
                         krange=nb.clust,
                         scaling=FALSE,
                         method=method.agnes)$stabk[nb.clust]}
  return(res.out)
}