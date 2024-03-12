#' @importFrom mix prelim.mix imp.mix em.mix
 
imputation.intern<-function(dataset,nb.clust,showits=TRUE, maxits=1000){
  flag<-FALSE
  if(is.data.frame(dataset)){
    flag<-TRUE
    don<-as.matrix(dataset)
  }else{don<-dataset}
  
  donclass<-cbind(class=c(nb.clust,rep(NA,nrow(don)-1)),don)
  # print(summary(donclass))
  s <- prelim.mix(donclass,1)
  # print(("ok"))
  thetahat <- try(em.mix(s,showits = showits,maxits = maxits),silent=TRUE)
  if(inherits(thetahat,"try-error")){warning(thetahat);return(NULL)}
  ximp <- imp.mix(s, thetahat, donclass)
  # print("ok3")
  clustering<-ximp[,1,drop=FALSE]
  ximp<-ximp[,-1]
  if(flag){ximp<-as.data.frame(ximp)}
  res.out<-list(ximp=ximp,cluster=clustering)
  return(res.out)
}