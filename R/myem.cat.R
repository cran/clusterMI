#' @importFrom cat prelim.cat em.cat

myem.cat<-function(dataset,K,boot=FALSE,silent=TRUE){
  # if(silent){sink("tmp.txt")}
  ind.boot<-seq(nrow(dataset))
  if(is.data.frame(dataset)){
    don<-as.matrix(dataset)
  }else{don<-dataset}
  don<-cbind(class=c(K,rep(NA,nrow(don)-1)),don)
  # print(don)
  if(!boot){s <- prelim.cat(don)
  }else{
    ind.boot<-sample(seq(nrow(don)),
                     size=nrow(don),
                     replace=TRUE)
    s <- prelim.cat(cbind(class=c(K,rep(NA,nrow(don)-1)),
                          don[ind.boot,]))
  }
  thetahat <- em.cat(s,showits=!silent)
  # if(silent){sink();file.remove("tmp.txt")}
  res.out<-list(s=s,theta=thetahat,ind.boot=ind.boot)
  return(res.out)
}