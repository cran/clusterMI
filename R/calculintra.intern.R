#' @importFrom FactoMineR MCA FAMD
#' @importFrom fpc kmeansCBI nselectboot claraCBI noisemclustCBI hclustCBI

calculintra.intern<-function(res.imp,method.clustering,Cboot,nb.clust,scaling,method.agnes,nstart.kmeans,modelNames,m.cmeans){
  
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
    # print(res.imp.intern)
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=clustermethod,
                         classification="centroid",krange=nb.clust,scaling=FALSE,count = FALSE,runs=nstart.kmeans)$stabk[nb.clust]}
  if(method.clustering=="cmeans"){
    clustermethod<-cmeansCBI.intern
    # print(res.imp.intern)
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=clustermethod,
                         classification="centroid",krange=nb.clust,scaling=FALSE,m.cmeans=m.cmeans)$stabk[nb.clust]}
  if(method.clustering%in%c("pam","clara")){
    clustermethod<-claraCBI
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=clustermethod,
                         classification="centroid",krange=nb.clust,usepam=(method.clustering=="pam"))$stabk[nb.clust]
  }
  if(method.clustering=="mixture"){
    # print(str(res.imp.intern))
    methodclassif<-ifelse(test=modelNames%in%c("EEE","EEI","EII"),yes = "lda",no="qda")
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=noisemclustCBI,
                         classification=methodclassif,krange=nb.clust,multipleboot=FALSE,modelNames=modelNames,verbose = FALSE)$stabk[nb.clust]
  }
  if(method.clustering=="agnes"){
    # print("okagnes")
    method.agnes.tmp<-method.agnes
    if(method.agnes=="gaverage"){method.agnes.tmp<-"average"}
    clustermethod<-hclustCBI
    res.out<-nselectboot(data=res.imp.intern,B=Cboot,clustermethod=clustermethod,
                         classification="averagedist",krange=nb.clust,scaling=FALSE,method=method.agnes.tmp)$stabk[nb.clust]}
  return(res.out)
}