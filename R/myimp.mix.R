#' @importFrom mix prelim.mix em.mix imp.mix rngseed
#' @import mix

myimp.mix<-function(dataset,K,seed=1234567,boot=FALSE,silent=TRUE,nbquali=0){
  # if(silent){sink("tmp.txt")}
  if(is.data.frame(dataset)){
    don<-as.matrix(dataset)
  }else{don<-dataset}
  
  # K<-3
  donclass<-cbind(class=c(K,rep(NA,nrow(don)-1)),don)
  s <- prelim.mix(donclass,1+nbquali)
  if(!boot){
    thetahat <- em.mix(s, showits = !silent)
  }else{
    ind.boot<-sample(seq(nrow(don)),
                     size=nrow(don),
                     replace=TRUE)
    donclassboot<-cbind(class=c(K,rep(NA,nrow(don)-1)),
                        don[ind.boot,])
    
    sboot <- prelim.mix(donclassboot,1+nbquali)
    # s$z<-donclass
    thetahat <- em.mix(sboot,showits = !silent)
  }
  
  
  if(!is.null(seed)){mix::rngseed(seed)}     # set random number generator seed
  ximp <- imp.mix(s, thetahat, donclass)
  # if(silent){sink();file.remove("tmp.txt")}
  if(!boot){
    return(ximp)
  }else{
    res.out<-list(ximp=ximp,theta=thetahat,s=s)
    return(res.out)
  }
}