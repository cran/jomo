jomo1con<- function(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, nburn=100, nbetween=100, nimp=5,output=1, out.iter=10) {
  if (nimp<2) {
    nimp=2
    cat("Minimum number of imputations:2. For single imputation using function jomo1con.MCMCchain\n")
  }
  if (is.null(X)) X=matrix(1,nrow(Y),1)
  if (is.null(beta.start)) beta.start=matrix(0,ncol(X),ncol(Y))
  if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
  if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
  for (i in 1:ncol(X)) {
    if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
  }
  stopifnot(nrow(Y)==nrow(X), nrow(beta.start)==ncol(X), ncol(beta.start)==ncol(Y),nrow(l1cov.start)==ncol(l1cov.start), nrow(l1cov.start)==ncol(Y), nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.prior)==nrow(l1cov.start))
  betait=matrix(0,nrow(beta.start),ncol(beta.start))
  for (i in 1:nrow(beta.start)) {
    for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
  }
  covit=matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  for (i in 1:nrow(l1cov.start)) {
    for (j in 1:ncol(l1cov.start)) covit[i,j]=l1cov.start[i,j]
  }   
  colnamy<-colnames(Y)
  colnamx<-colnames(X)
  Y<-data.matrix(Y)
  storage.mode(Y) <- "numeric"    
  X<-data.matrix(X)
  storage.mode(X) <- "numeric"    
  if (output!=1) out.iter=nburn+nbetween
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+2)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(X), (ncol(Y)+ncol(X)+1)]=c(1:nrow(Y))
  Yimp=Y
  Yimp2=matrix(Yimp, nrow(Y),ncol(Y))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+1)]=c(1:nrow(Y))
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+2)]=1
  betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),(nimp-1)))
  bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
  omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),(nimp-1)))
  opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  meanobs<-colMeans(Y,na.rm=TRUE)
  for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
  Y.cat<-Y.numcat<-(-999)
  .Call("jomo1C", Y, Yimp, Yimp2, Y.cat, X,betait,bpost,covit,opost, nburn, l1cov.prior,Y.numcat, ncol(Y),out.iter,0, PACKAGE = "jomo")
  #betapost[,,1]=bpost
  #omegapost[,,1]=opost
  bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
  opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
  if (output==1) cat("First imputation registered.", "\n")
  for (i in 2:nimp) {
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+1)]=c(1:nrow(Y))
    imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+2)]=i
    .Call("jomo1C", Y, Yimp, Yimp2, Y.cat, X,betait,bpost,covit, opost, nbetween, l1cov.prior, Y.numcat, ncol(Y),out.iter,0, PACKAGE = "jomo")
    betapost[,,(i-1)]=bpost
    omegapost[,,(i-1)]=opost
    bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
    opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    imp[(i*nrow(Y)+1):((i+1)*nrow(Y)),1:ncol(Y)]=Yimp2
    if (output==1) cat("Imputation number ", i, "registered", "\n")
  }
  imp<-data.frame(imp)
  if (is.null(colnamy)) colnamy=paste("Y", 1:ncol(Y), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
  dimnames(betapost)[1] <- list(colnamx)
  dimnames(betapost)[2] <- list(colnamy)
  dimnames(omegapost)[1] <- list(colnamy)
  dimnames(omegapost)[2] <- list(colnamy)
  betapostmean<-apply(betapost, c(1,2), mean)
  omegapostmean<-apply(omegapost, c(1,2), mean)
  if (output==1) {
    cat("The posterior mean of the fixed effects estimates is:\n")
    print(t(betapostmean))
    cat("\nThe posterior covariance matrix is:\n")
    print(omegapostmean)
  }
  colnames(imp)<-c(colnamy,colnamx,"id","Imputation")
  return(imp)
}