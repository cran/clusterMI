drawW<-function(res.mclust,method="QDA"){
  # save(res.mclust,file="res.mclust.Rdata")
  # print(summary(res.mclust$data))#;print((res.mclust$G));print((res.mclust$parameters))
  #res.mclust contient parfois des sorties incohérente (les parametres ne correspondent pas au nb de groupes !)
  if(length(res.mclust$parameters$pro)!=res.mclust$G){
    warning("Mclust: the number of components != length of proportion. Model is fitted again")
    res.mclust.tmp<-mclust::Mclust(res.mclust$data,G=res.mclust$G,modelNames=res.mclust$modelName,verbose = FALSE)
    if(length(res.mclust.tmp$parameters$pro)!=res.mclust.tmp$G){stop("issue in Mclust: the number of components != length of proportion")}
  }else{res.mclust.tmp<-res.mclust}
  nb.clust<-res.mclust.tmp$G
  n<-length(res.mclust.tmp$classification)
  mu<-lapply(as.data.frame(res.mclust.tmp$parameters$mean),I)
  theta<-as.list(res.mclust.tmp$parameters$pro)
  Sigma<-list();for(ii in seq(nb.clust)){Sigma[[ii]]<-res.mclust.tmp$parameters$variance$sigma[,,ii]}
  Z<-res.mclust.tmp$data
  if(method=="QDA"){
    proba<- mapply(FUN = function(mu,theta,Sigma,Z){
      Zcentre<-Z-matrix(mu,nrow(Z),ncol(Z),byrow=TRUE)
      cst<- -.5*log(det(Sigma))+log(theta)
      res.out<-exp(-0.5*diag(tcrossprod(Zcentre%*%solve(Sigma),Zcentre))+cst)
      
      return(res.out)
    },mu=mu,theta=theta,Sigma=Sigma,MoreArgs = list(Z=Z))
  }else if(method=="LDA"){
    proba<- mapply(FUN = function(mu,theta,Sigma,Z){
      sigmainvmu<-solve(Sigma)%*%matrix(mu,ncol=1)
      cst<- (-.5*(matrix(mu,nrow=1)%*%sigmainvmu)+log(theta))
      res.out<-t(exp(Z%*%sigmainvmu+rep(cst[1,1],nrow(Z))))
      return(res.out)
    },mu=mu,theta=theta,Sigma=Sigma,MoreArgs = list(Z=Z))
    # print(proba)
  }
  
  W<-apply(proba,1,FUN = function(xx){
    if(any(is.infinite(xx))){
      #tirage au hasard si proba= Inf
      res.out<-sample(seq(length(xx)),size=1)
      warning("Inf in probabilities")
    }else{
      
      res.out<-sample(seq(length(xx)),prob = xx,size=1)}
    return(res.out)
  })
  return(W)
}