jomo1ranmixhr <-
function(Y.con, Y.cat, Y.numcat, X=NULL, Z=NULL, clus, beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, nburn=1000, nbetween=1000, nimp=5, a=NULL, a.prior=NULL, meth="random", output=1, out.iter=10) {
  if (nimp<2) {
    nimp=2
    cat("Minimum number of imputations:2. For single imputation using function jomo1ranmixhr.MCMCchain\n")
  }
  if (is.null(X)) X=matrix(1,nrow(Y.cat),1) 
  if (is.null(Z)) Z=matrix(1,nrow(Y.cat),1)
  if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
  if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
  if (is.null(a)) a=ncol(beta.start)+50
  if (is.null(a.prior)) a.prior=ncol(beta.start)
  clus<-factor(unlist(clus))
  previous_levels_clus<-levels(clus)
  levels(clus)<-0:(nlevels(clus)-1)
  if (is.null(u.start)) u.start = matrix(0, nlevels(clus),ncol(Z)*(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
  if (is.null(l2cov.start)) l2cov.start = diag(1, ncol(u.start))
  if (is.null(l2cov.prior)) l2cov.prior = diag(1, ncol(l2cov.start))
  if (is.null(l1cov.start)) l1cov.start=matrix(diag(1,ncol(beta.start)),ncol(beta.start)*nlevels(clus),ncol(beta.start),2)
  previous_levels<-list()
  Y.cat<-data.frame(Y.cat)
  for (i in 1:ncol(Y.cat)) {
    Y.cat[,i]<-factor(Y.cat[,i])
    previous_levels[[i]]<-levels(Y.cat[,i])
    levels(Y.cat[,i])<-1:nlevels(Y.cat[,i])
  }
  for (i in 1:ncol(X)) {
    if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
  }
  for (i in 1:ncol(Z)) {
    if (is.factor(Z[,i])) Z[,i]<-as.numeric(Z[,i])
  }
  stopifnot((meth=="fixed"|meth=="random"),nrow(Y.con)==nrow(clus),nrow(Y.con)==nrow(X), nrow(beta.start)==ncol(X), ncol(beta.start)==(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))),nrow(l1cov.start)==nrow(u.start)*ncol(l1cov.start), nrow(l1cov.start)==nrow(u.start)*ncol(beta.start), nrow(l1cov.prior)==ncol(l1cov.prior),nrow(l1cov.start)==nrow(u.start)*nrow(l1cov.prior),nrow(Z)==nrow(Y.con), ncol(l2cov.start)==ncol(u.start), ncol(u.start)==ncol(Z)*(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
  betait=matrix(0,nrow(beta.start),ncol(beta.start))
  for (i in 1:nrow(beta.start)) {
    for (j in 1:ncol(beta.start)) betait[i,j]=beta.start[i,j]
  }
  covit=matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  for (i in 1:nrow(l1cov.start)) {
    for (j in 1:ncol(l1cov.start)) covit[i,j]=l1cov.start[i,j]
  }   
  uit=matrix(0,nrow(u.start),ncol(u.start))
  for (i in 1:nrow(u.start)) {
    for (j in 1:ncol(u.start)) uit[i,j]=u.start[i,j]
  }
  covuit=matrix(0,nrow(l2cov.start),ncol(l2cov.start))
  for (i in 1:nrow(l2cov.start)) {
    for (j in 1:ncol(l2cov.start)) covuit[i,j]=l2cov.start[i,j]
  }   
  ait=as.numeric(a)
  colnamycon<-colnames(Y.con)
  colnamycat<-colnames(Y.cat)
  colnamx<-colnames(X)
  colnamz<-colnames(Z)
  Y.con<-data.matrix(Y.con)
  storage.mode(Y.con) <- "numeric"    
  Y.cat<-data.matrix(Y.cat)
  storage.mode(Y.cat) <- "numeric"    
  X<-data.matrix(X)
  storage.mode(X) <- "numeric"  
  stopifnot(!any(is.na(X)))
  Z<-data.matrix(Z)
  storage.mode(Z) <- "numeric"
  stopifnot(!any(is.na(Z)))
  clus <- matrix(as.integer(levels(clus))[clus], ncol=1)
  Y=cbind(Y.con,Y.cat)
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
  Yi=cbind(Y.con, matrix(0,nrow(Y.con),(sum(Y.numcat)-length(Y.numcat))))
  h=1
  for (i in 1:length(Y.numcat)) {
    for (j in 1:nrow(Y)) {
      if (is.na(Y.cat[j,i])) {
        Yi[j,(ncol(Y.con)+h):(ncol(Y.con)+h+Y.numcat[i]-2)]=NA
      }
    } 
    h=h+Y.numcat[i]-1
  }
  if (output!=1) out.iter=nburn+nbetween
  imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+ncol(Z)+3)
  imp[1:nrow(Y),1:ncol(Y)]=Y
  imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[1:nrow(Z), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[1:nrow(clus), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[1:nrow(X), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  Yimp=Yi
  Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
  imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
  imp[(nrow(Z)+1):(2*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
  imp[(nrow(clus)+1):(2*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
  imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=1  
  betapost<- array(0, dim=c(nrow(beta.start),ncol(beta.start),(nimp-1)))
  bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
  upost<-matrix(0,nrow(u.start),ncol(u.start))
  upostall<-array(0, dim=c(nrow(u.start),ncol(u.start),(nimp-1)))
  omegapost<- array(0, dim=c(nrow(l1cov.start),ncol(l1cov.start),(nimp-1)))
  opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  covupost<- array(0, dim=c(nrow(l2cov.start),ncol(l2cov.start),(nimp-1)))
  cpost<-matrix(0,nrow(l2cov.start),ncol(l2cov.start))
  meanobs<-colMeans(Yi,na.rm=TRUE)
  for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=rnorm(1,meanobs[j],1)
  if (meth=="fixed") {
    fixed=1    
  } else {
    fixed=0
  }
    .Call("jomo1ranhrC", Y, Yimp, Yimp2, Y.cat, X, Z, clus,betait,uit,bpost,upost,covit,opost, covuit,cpost,nburn, l1cov.prior,l2cov.prior,Y.numcat, ncol(Y.con),ait, a.prior, out.iter, fixed, 0, miss.pat.id, n.patterns, PACKAGE = "jomo")
  #betapost[,,1]=bpost
  #upostall[,,1]=upost
  #omegapost[,,(1)]=opost
  #covupost[,,(1)]=cpost  
  bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
  opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
  upost<-matrix(0,nrow(u.start),ncol(u.start))
  cpost<-matrix(0,nrow(l2cov.start),ncol(l2cov.start))
  imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y.con)]=Yimp2[,1:ncol(Y.con)]
  imp[(nrow(Y)+1):(2*nrow(Y)),(ncol(Y.con)+1):ncol(Y)]=Y.cat
  if (output==1) cat("First imputation registered.", "\n")
  for (i in 2:nimp) {
    #Yimp2=matrix(0, nrow(Yimp),ncol(Yimp))
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
    imp[(i*nrow(clus)+1):((i+1)*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
    imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=i  
    if (meth=="fixed") {
      fixed=1    
    } else {
      fixed=0
    }
    .Call("jomo1ranhrC", Y, Yimp, Yimp2, Y.cat, X, Z, clus,betait,uit,bpost,upost,covit,opost, covuit,cpost,nbetween, l1cov.prior,l2cov.prior,Y.numcat, ncol(Y.con),ait,a.prior,out.iter, fixed, 0, miss.pat.id, n.patterns, PACKAGE = "jomo")
    betapost[,,(i-1)]=bpost
    upostall[,,(i-1)]=upost
    omegapost[,,(i-1)]=opost
    covupost[,,(i-1)]=cpost
    bpost<-matrix(0,nrow(beta.start),ncol(beta.start))
    opost<-matrix(0,nrow(l1cov.start),ncol(l1cov.start))
    upost<-matrix(0,nrow(u.start),ncol(u.start))
    cpost<-matrix(0,nrow(l2cov.start),ncol(l2cov.start))
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),1:ncol(Y.con)]=Yimp2[,1:ncol(Y.con)]
    imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y.con)+1):ncol(Y)]=Y.cat
    if (output==1) cat("Imputation number ", i, "registered", "\n")
  }
  
  imp<-data.frame(imp)
  for (i in 1:ncol(Y.cat)) {
    imp[,(ncol(Y.con)+i)]<-as.factor(imp[,(ncol(Y.con)+i)]) 
    levels(imp[,(ncol(Y.con)+i)])<-previous_levels[[i]]
  }
  imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)]<-factor(imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)])
  levels(imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)])<-previous_levels_clus
  clus<-factor(clus)
  levels(clus)<-previous_levels_clus
  for (j in 1:(ncol(Y.con))) {
    imp[,j]=as.numeric(imp[,j])
  }
  for (j in (ncol(Y)+1):(ncol(Y)+ncol(X)+ncol(Z))) {
    imp[,j]=as.numeric(imp[,j])
  }
  if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y.cat), sep = "")
  if (is.null(colnamycon)) colnamycon=paste("Ycon", 1:ncol(Y.con), sep = "")
  if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
  if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")  
  colnames(imp)<-c(colnamycon,colnamycat,colnamx,colnamz,"clus","id","Imputation")
  cnycatcomp<-rep(NA,(sum(Y.numcat)-length(Y.numcat)))
  count=0
  for ( j in 1:ncol(Y.cat)) {
    for (k in 1:(Y.numcat[j]-1)) {
      cnycatcomp[count+k]<-paste(colnamycat[j],k,sep=".")
    }
    count=count+Y.numcat[j]-1
  }
  cnamycomp<-c(colnamycon,cnycatcomp)
  dimnames(betapost)[1] <- list(colnamx)
  dimnames(betapost)[2] <- list(cnamycomp)
  dimnames(omegapost)[1] <- list(paste(cnamycomp,rep(levels(clus),each=ncol(Yimp2)), sep="."))
  dimnames(omegapost)[2] <- list(cnamycomp)
  colnamcovu<-paste(cnamycomp,rep(colnamz,each=ncol(omegapost)),sep="*")
  dimnames(covupost)[1] <- list(colnamcovu)
  dimnames(covupost)[2] <- list(colnamcovu)
  dimnames(upostall)[1]<-list(levels(clus))
  dimnames(upostall)[2]<-list(colnamcovu)
  betapostmean<-data.frame(apply(betapost, c(1,2), mean))    
  upostmean<-data.frame(apply(upostall, c(1,2), mean))
  omegapostmean<-data.frame(apply(omegapost, c(1,2), mean))
  covupostmean<-data.frame(apply(covupost, c(1,2), mean))
  if (output==1) {
    cat("The posterior mean of the fixed effects estimates is:\n")
    print(t(betapostmean))
    cat("\nThe posterior mean of the random effects estimates is:\n")
    print(upostmean)
    cat("\nThe posterior mean of the level 1 covariance matrices is:\n")
    print(omegapostmean)
    cat("\nThe posterior mean of the level 2 covariance matrix is:\n")
    print(covupostmean)
  }
  return(imp)
}
