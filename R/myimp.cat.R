#' @importFrom cat prelim.cat em.cat imp.cat

myimp.cat<-function(dataset,K,seed=1234567,boot=FALSE,silent=TRUE){
  # if(silent){sink("tmp.txt")}
  if(is.data.frame(dataset)){
    don<-as.matrix(dataset)
  }else{don<-dataset}
  
  # K<-3
  donclass<-cbind(class=c(K,rep(NA,nrow(don)-1)),don)
  s <- prelim.cat(donclass)
  if(!boot){
    thetahat <- em.cat(s)
  }else{
    ind.boot<-sample(seq(nrow(don)),
                     size=nrow(don),
                     replace=TRUE)
    donclassboot<-cbind(class=c(K,rep(NA,nrow(don)-1)),
                        don[ind.boot,])
    
    sboot <- prelim.cat(donclassboot)
    # s$z<-donclass
    thetahat <- em.cat(sboot,showits = !silent)
  }
  
  
  if(!is.null(seed)){rngseedcat(seed)}     # set random number generator seed
  ximp <- imp.cat(s, thetahat)
  colnames(ximp)<-c("class",colnames(dataset))
  # if(silent){sink();file.remove("tmp.txt")}
  if(!boot){
    return(ximp)
  }else{
    res.out<-list(ximp=ximp,theta=thetahat,s=s)
    return(res.out)
  }
}