#' @importFrom mix prelim.mix em.mix

myem.mix<-function(dataset,K,boot=FALSE,silent=TRUE,nbquali=0){
  # if(silent){sink("tmp.txt")}
  ind.boot<-seq(nrow(dataset))
  if(is.data.frame(dataset)){
    don<-as.matrix(dataset)
  }else{don<-dataset}
  don<-cbind(class=c(K,rep(NA,nrow(don)-1)),don)
  # print(don)
  if(!boot){s <- prelim.mix(don,1+nbquali)
  }else{
    ind.boot<-sample(seq(nrow(don)),
                     size=nrow(don),
                     replace=TRUE)
    s <- prelim.mix(cbind(class=c(K,rep(NA,nrow(don)-1)),
                               don[ind.boot,]),1+nbquali)
  }
  thetahat <- em.mix(s,showits=!silent)
  # if(silent){sink();file.remove("tmp.txt")}
  res.out<-list(s=s,theta=thetahat,ind.boot=ind.boot)
  return(res.out)
}