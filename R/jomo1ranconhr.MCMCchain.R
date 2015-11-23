jomo1ranconhr.MCMCchain <-
  function(Y, X=matrix(1,nrow(Y),1), Z=matrix(1,nrow(Y),1), clus, betap=matrix(0,ncol(X),ncol(Y)), up=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y)), covp=matrix(diag(1,ncol(Y)),nrow(unique(clus))*ncol(Y),ncol(Y),2), covu=diag(1,ncol(Y)*ncol(Z)), Sp=diag(1,ncol(Y)), Sup=diag(1,ncol(Y)*ncol(Z)), nburn=100, a=ncol(Y),meth="random", output=1, out.iter=10) {
    clus<-factor(unlist(clus))
    previous_levels_clus<-levels(clus)
    levels(clus)<-0:(nlevels(clus)-1)
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    for (i in 1:ncol(Z)) {
      if (is.factor(Z[,i])) Z[,i]<-as.numeric(Z[,i])
    }
    stopifnot((meth=="fixed"|meth=="random"),nrow(Y)==nrow(clus),nrow(Y)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==ncol(Y),nrow(covp)==nrow(up)*ncol(covp), nrow(covp)==nrow(up)*ncol(Y), nrow(Sp)==ncol(Sp),nrow(covp)==nrow(up)*nrow(Sp), nrow(Z)==nrow(Y), ncol(covu)==ncol(up), ncol(up)==ncol(Z)*ncol(Y))
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
    ait=0
    ait=a
    nimp=1
    colnamy<-colnames(Y)
    colnamx<-colnames(X)
    colnamz<-colnames(Z)
    Y<-data.matrix(Y)
    storage.mode(Y) <- "numeric"    
    X<-data.matrix(X)
    storage.mode(X) <- "numeric"    
    Z<-data.matrix(Z)
    storage.mode(Z) <- "numeric"    
    clus<-data.matrix(clus)
    if (output!=1) out.iter=nburn+2
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
    upostall<-array(0, dim=c(nrow(up),ncol(up),nburn))
    covupost<- array(0, dim=c(nrow(covu),ncol(covu),nburn))
    meanobs<-colMeans(Y,na.rm=TRUE)
    for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=meanobs[j]
    #for (i in 1:nrow(Y)) for (j in 1:ncol(Y)) if (is.na(Yimp[i,j])) Yimp[i,j]=rnorm(1,mean=meanobs[j], sd=0.01)
    if (meth=="fixed") {
      .Call("MCMCjomo1ranconhf", Y, Yimp, Yimp2, X, Z, clus, betait, uit, betapost, upostall, covit,omegapost, covuit,covupost, nburn, Sp, Sup,out.iter,  PACKAGE = "jomo") 
    }
    if (meth=="random") {
      .Call("MCMCjomo1ranconhr", Y, Yimp, Yimp2, X, Z, clus, betait, uit, betapost, upostall, covit,omegapost, covuit,covupost, nburn, Sp, Sup, ait,out.iter, PACKAGE = "jomo") 
    }
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Yimp2
    Yimp=Yimp2
    betapostmean<-apply(betapost, c(1,2), mean)
    upostmean<-apply(upostall, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    covupostmean<-apply(covupost, c(1,2), mean)
    if (output==1) {
      cat("The posterior mean of the fixed effects estimates is:\n")
      print(betapostmean)
      cat("The posterior mean of the random effects estimates is:\n")
      print(upostmean)
      cat("The posterior mean of the level 1 covariance matrices is:\n")
      print(omegapostmean)
      cat("The posterior mean of the level 2 covariance matrix is:\n")
      print(covupostmean)
    }
    imp<-data.frame(imp)
    levels(imp[,(ncol(Y)+ncol(X)+ncol(Z)+1)])<-previous_levels_clus
    levels(clus)<-previous_levels_clus
    if (is.null(colnamy)) colnamy=paste("Y", 1:ncol(Y), sep = "")
    if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamy,colnamx,colnamz,"clus","id","Imputation")
    return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost,"collectu"=upostall, "collectcovu"=covupost))
  }
