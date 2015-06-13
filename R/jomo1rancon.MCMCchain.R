jomo1rancon.MCMCchain<- function(Y, X=matrix(1,nrow(Y),1), Z=matrix(1,nrow(Y),1), clus, betap=matrix(0,ncol(X),ncol(Y)), up=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y)), covp=diag(1,ncol(Y)), covu=diag(1,ncol(Y)*ncol(Z)), Sp=diag(1,ncol(Y)), Sup=diag(1,ncol(Y)*ncol(Z)), nburn=100) {
  stopifnot(nrow(Y)==nrow(clus),nrow(Y)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==ncol(Y),nrow(covp)==ncol(covp), nrow(covp)==ncol(Y), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp), nrow(Z)==nrow(Y), ncol(covu)==ncol(up), ncol(up)==ncol(Z)*ncol(Y))
  betait=matrix(0,nrow(betap),ncol(betap))
  for (i in 1:nrow(betap)) {
    for (j in 1:ncol(betap)) betait[i,j]=betap[i,j]
  }
  covit=matrix(0,nrow(covp),ncol(covp))
  for (i in 1:nrow(covp)) {
    for (j in 1:ncol(covp)) covit[i,j]=covp[i,j]
  }   
  uit=matrix(0,nrow(up),ncol(up))
  for (i in 1:nrow(up)) {
    for (j in 1:ncol(up)) uit[i,j]=up[i,j]
  }
  covuit=matrix(0,nrow(covu),ncol(covu))
  for (i in 1:nrow(covu)) {
    for (j in 1:ncol(covu)) covuit[i,j]=covu[i,j]
  }   
  rngflag=0
  nimp=1
  colnamy<-colnames(Y)
  colnamx<-colnames(X)
  colnamz<-colnames(Z)
  Y<-as.matrix(Y,nrow(Y),ncol(Y))
  X<-as.matrix(X,nrow(X),ncol(X))
  Z<-as.matrix(Z,nrow(Z),ncol(Z))
  clus<-as.matrix(clus,nrow(clus),ncol(clus))
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+ncol(Z)+3)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(Z), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[1:nrow(clus), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[1:nrow(X), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  Yimp=Y
  Yimp2=matrix(0, nrow(Y),ncol(Y))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(Z)+1):(2*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[(nrow(clus)+1):(2*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=1
  betapost<- array(0, dim=c(nrow(betap),ncol(betap),nburn))
  omegapost<- array(0, dim=c(nrow(covp),ncol(covp),nburn))
  collectimp<- array(0, dim=c(nrow(Y),ncol(covp),nburn))
  upostall<-array(0, dim=c(nrow(up),ncol(up),nburn))
  covupost<- array(0, dim=c(nrow(covu),ncol(covu),nburn))
  meanobs<-colMeans(Y,na.rm=TRUE)
  varobs<-apply(Y,2,sd,na.rm=TRUE)
  for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=meanobs[j]
  #for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=rnorm(1,mean=meanobs[j], sd=0.01)
  .Call("MCMCjomo1rancon", Y, Yimp, Yimp2, X, Z, clus, betait, uit, betapost, upostall, covit, omegapost, covuit, covupost, nburn, Sp, Sup,rngflag,collectimp, PACKAGE = "jomo")
 
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
  Yimp=Yimp2
  betapostmean<-apply(betapost, c(1,2), mean)
  upostmean<-apply(upostall, c(1,2), mean)
  omegapostmean<-apply(omegapost, c(1,2), mean)
  covupostmean<-apply(covupost, c(1,2), mean)
  cat("The posterior mean of the fixed effects estimates is:\n")
  print(betapostmean)
  cat("The posterior mean of the random effects estimates is:\n")
  print(upostmean)
  cat("The posterior mean of the level 1 covariance matrix is:\n")
  print(omegapostmean)
  cat("The posterior mean of the level 2 covariance matrix is:\n")
  print(covupostmean)
  imp<-data.frame(imp)
  if (is.null(colnamy)) colnamy=paste("Y", 1:ncol(Y), sep = "")
  if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
  colnames(imp)<-c(colnamy,colnamx,colnamz,"clus","id","Imputation")
  return(list("finimp"=imp,"collectimp"=collectimp,"collectbeta"=betapost,"collectomega"=omegapost,"collectu"=upostall, "collectcovu"=covupost))
}