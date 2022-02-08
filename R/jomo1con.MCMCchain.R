jomo1con.MCMCchain<- function(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL,start.imp=NULL, nburn=100,output=1, out.iter=10) {
  if (is.null(X)) X=matrix(1,nrow(Y),1)
  if (is.null(beta.start)) beta.start=matrix(0,ncol(X),ncol(Y))
  if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
  if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
  if (is_tibble(Y)) {
    Y<-data.frame(Y)
    warning("tibbles not supported. Y converted to standard data.frame. ")
  }
  if (is_tibble(X)) {
    X<-data.frame(X)
    warning("tibbles not supported. X converted to standard data.frame. ")
  }
  
  if (any(is.na(Y))) {
    if (ncol(Y)==1) {
      miss.pat<-matrix(c(0,1),2,1)
      n.patterns<-2
    } else  {
      miss.pat<-md.pattern.mice(Y, plot=F)
      miss.pat<-miss.pat[,colnames(Y)]
      n.patterns<-nrow(miss.pat)-1
    }
  } else {
    miss.pat<-matrix(0,2,ncol(Y)+1)
    n.patterns<-nrow(miss.pat)-1
  }
  
  miss.pat.id<-rep(0,nrow(Y))
  for (i in 1:nrow(Y)) {
    k <- 1
    flag <- 0
    while ((k <= n.patterns) & (flag == 0)) {
      if (all(!is.na(Y[i,])==miss.pat[k,1:(ncol(miss.pat))])) {
        miss.pat.id[i] <- k
        flag <- 1
      } else {
        k <- k + 1
      }
    }
  }
  
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
  nimp=1
  colnamy<-colnames(Y)
  colnamx<-colnames(X)
  Y<-data.matrix(Y)
  storage.mode(Y) <- "numeric"    
  X<-data.matrix(X)
  storage.mode(X) <- "numeric"  
  stopifnot(!any(is.na(X)))
  if (output!=1) out.iter=nburn+2
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+2)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(X), (ncol(Y)+ncol(X)+1)]=c(1:nrow(Y))
  Yimp=Y
  Yimp2=matrix(Yimp, nrow(Y),ncol(Y))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+1)]=c(1:nrow(Y))
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+2)]=1
  betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),nburn))
  omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),nburn))
  meanobs<-colMeans(Y,na.rm=TRUE)
  if (!is.null(start.imp)) {
    start.imp<-as.matrix(start.imp)
    if ((nrow(start.imp)!=nrow(Yimp2))||(ncol(Yimp2)>ncol(start.imp))) {
      cat("start.imp dimensions incorrect. Not using start.imp as starting value for the imputed dataset.\n")
      start.imp=NULL
    } else {
      if ((nrow(start.imp)==nrow(Yimp2))&(ncol(Yimp2)<ncol(start.imp))) {
        Yimp2<-start.imp[,1:ncol(Yimp2)]
        cat("NOTE: start.imp has more columns than needed. Dropping unnecessary columns.\n")
      } else {
        Yimp2<-start.imp
      }
    }
  }
  if (is.null(start.imp)) {
    for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
  } 
  Y.cat<-Y.numcat<-(-999)
  .Call("jomo1C", Y, Yimp, Yimp2, Y.cat, X,betait,betapost,covit,omegapost, nburn, l1cov.prior,Y.numcat, ncol(Y),out.iter, 1, miss.pat.id, n.patterns, PACKAGE = "jomo")
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
  imp<-data.frame(imp)
  if (is.null(colnamy)) colnamy=paste("Y", 1:ncol(Y), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
  colnames(imp)<-c(colnamy,colnamx,"id","Imputation")
  dimnames(betapost)[1] <- list(colnamx)
  dimnames(betapost)[2] <- list(colnamy)
  dimnames(omegapost)[1] <- list(colnamy)
  dimnames(omegapost)[2] <- list(colnamy)
  betapostmean<-data.frame(apply(betapost, c(1,2), mean))
  omegapostmean<-data.frame(apply(omegapost, c(1,2), mean))
  if (output==1) {
    cat("The posterior mean of the fixed effects estimates is:\n")
    print(t(betapostmean))
    cat("\nThe posterior covariance matrix is:\n")
    print(omegapostmean)
  }
  return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost))
}