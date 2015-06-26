jomo1con.MCMCchain<- function(Y, X=matrix(1,nrow(Y),1), betap=matrix(0,ncol(X),ncol(Y)), covp=diag(1,ncol(Y)), Sp=diag(1,ncol(Y)), nburn=100,output=1, out.iter=10) {
  for (i in 1:ncol(X)) {
    if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
  }
  stopifnot(nrow(Y)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==ncol(Y),nrow(covp)==ncol(covp), nrow(covp)==ncol(Y), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp))
  betait=matrix(0,nrow(betap),ncol(betap))
  for (i in 1:nrow(betap)) {
    for (j in 1:ncol(betap)) betait[i,j]=betap[i,j]
  }
  covit=matrix(0,nrow(covp),ncol(covp))
  for (i in 1:nrow(covp)) {
    for (j in 1:ncol(covp)) covit[i,j]=covp[i,j]
  }   
  nimp=1
  colnamy<-colnames(Y)
  colnamx<-colnames(X)
  Y<-as.matrix(Y,nrow(Y),ncol(Y))
  X<-as.matrix(X,nrow(X),ncol(X))
  if (output!=1) out.iter=nburn+2
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+2)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(X), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
  Yimp=Y
  Yimp2=matrix(0, nrow(Y),ncol(Y))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+1)]=1
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
  betapost<- array(0, dim=c(nrow(betap),ncol(betap),nburn))
  omegapost<- array(0, dim=c(nrow(covp),ncol(covp),nburn))
  meanobs<-colMeans(Y,na.rm=TRUE)
  for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=meanobs[j]
  .Call("MCMCjomo1con", Y, Yimp, Yimp2, X,betait,betapost,covit, omegapost, nburn, Sp,out.iter, PACKAGE = "jomo")
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
  Yimp=Yimp2
  betapostmean<-apply(betapost, c(1,2), mean)
  omegapostmean<-apply(omegapost, c(1,2), mean)
  if (output==1) {
    cat("The posterior mean of the fixed effects estimates is:\n")
    print(betapostmean)
    cat("The posterior covariance matrix is:\n")
    print(omegapostmean)
  }
  imp<-data.frame(imp)
  if (is.null(colnamy)) colnamy=paste("Y", 1:ncol(Y), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
  colnames(imp)<-c(colnamy,colnamx,"Imputation","id")
  return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost))
}