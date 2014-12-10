jomo1rancon<- function(Y, X=matrix(1,nrow(Y),1), Z=matrix(1,nrow(Y),1), clus, betap=matrix(0,ncol(X),ncol(Y)), up=matrix(0,length(unique(clus)),ncol(Z)*ncol(Y)), covp=diag(1,ncol(Y)), covu=diag(1,ncol(Y)*ncol(Z)), Sp=diag(1,ncol(Y)), Sup=diag(1,ncol(Y)*ncol(Z)), nburn=100, nbetween=100, nimp=5, rngflag=0) {
  stopifnot(nrow(Y)==nrow(clus),nrow(Y)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==ncol(Y),nrow(covp)==ncol(covp), nrow(covp)==ncol(Y), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp), nrow(Z)==nrow(Y), ncol(covu)==ncol(up), ncol(up)==ncol(Z)*ncol(Y))
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
  betapost<- array(0, dim=c(nrow(betap),ncol(betap),(nimp-1)))
  bpost<-matrix(0,nrow(betap),ncol(betap))
  upost<-matrix(0,nrow(up),ncol(up))
  upostall<-array(0, dim=c(nrow(up),ncol(up),(nimp-1)))
  omegapost<- array(0, dim=c(nrow(covp),ncol(covp),(nimp-1)))
  opost<-matrix(0,nrow(covp),ncol(covp))
  covupost<- array(0, dim=c(nrow(covu),ncol(covu),(nimp-1)))
  cpost<-matrix(0,nrow(covu),ncol(covu))
  meanobs<-colMeans(Y,na.rm=TRUE)
  varobs<-apply(Y,2,sd,na.rm=TRUE)
  for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=meanobs[j]
  #for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=rnorm(1,mean=meanobs[j], sd=0.01)
  .Call("jomo1rancon", Y, Yimp, Yimp2, X, Z, clus, betap, up, bpost, upost, covp, opost, covu, cpost, nburn, Sp, Sup,rngflag, PACKAGE = "jomo")
  #betapost[,,1]=bpost
  #upostall[,,1]=upost
  #omegapost[,,1]=opost
  #covupost[,,1]=cpost
  bpost<-matrix(0,nrow(betap),ncol(betap))
  opost<-matrix(0,nrow(covp),ncol(covp))
  upost<-matrix(0,nrow(up),ncol(up))
  cpost<-matrix(0,nrow(covu),ncol(covu))
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
  Yimp=Yimp2
  cat("First imputation registered.", "\n")
  for (i in 2:nimp) {
    Yimp2=matrix(0, nrow(Y),ncol(Y))
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
    imp[(i*nrow(clus)+1):((i+1)*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=i
    .Call("jomo1rancon", Y, Yimp, Yimp2, X, Z, clus, betap, up, bpost, upost, covp,opost, covu,cpost, nbetween, Sp, Sup,rngflag, PACKAGE = "jomo")
    betapost[,,(i-1)]=bpost
    upostall[,,(i-1)]=upost
    omegapost[,,(i-1)]=opost
    covupost[,,(i-1)]=cpost
    bpost<-matrix(0,nrow(betap),ncol(betap))
    opost<-matrix(0,nrow(covp),ncol(covp))
    upost<-matrix(0,nrow(up),ncol(up))
    cpost<-matrix(0,nrow(covu),ncol(covu))
    imp[(i*nrow(Y)+1):((i+1)*nrow(Y)),1:ncol(Y)]=Yimp2
    Yimp=Yimp2
    cat("Imputation number ", i, "registered", "\n")
  }
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
  colnames(imp)<-c(colnamy,colnamx,colnamz,"clus","id","Imputation")
  return(imp)
}