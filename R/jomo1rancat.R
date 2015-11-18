jomo1rancat <-
  function(Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), Z=matrix(1,nrow(Y_cat),1), clus, betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat)))), up=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y_numcat)-length(Y_numcat)))), covp=diag(1,ncol(betap)), covu=diag(1,ncol(up)), Sp=diag(1,ncol(covp)), Sup=diag(1,ncol(covu)), nburn=100, nbetween=100, nimp=5, output=1, out.iter=10) {
    previous_levels<-list()
    for (i in 1:ncol(Y_cat)) {
      Y_cat[,i]<-factor(Y_cat[,i])
      previous_levels[[i]]<-levels(Y_cat[,i])
      levels(Y_cat[,i])<-1:nlevels(Y_cat[,i])
    }
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    for (i in 1:ncol(Z)) {
      if (is.factor(Z[,i])) Z[,i]<-as.numeric(Z[,i])
    }
    stopifnot( nrow(betap)==ncol(X), ncol(betap)==((sum(Y_numcat)-length(Y_numcat))),nrow(covp)==ncol(covp), nrow(covp)==ncol(betap), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp),nrow(Z)==nrow(Y_cat), ncol(covu)==ncol(up), ncol(up)==ncol(Z)*((sum(Y_numcat)-length(Y_numcat))))
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
    colnamycat<-colnames(Y_cat)
    colnamx<-colnames(X)
    colnamz<-colnames(Z)
    Y_cat<-data.matrix(Y_cat)
    X<-data.matrix(X)
    Z<-data.matrix(Z)
    clus<-data.matrix(clus)
    Y=cbind(Y_cat)
    Yi=cbind(matrix(0,nrow(Y_cat),(sum(Y_numcat)-length(Y_numcat))))
    h=1
    for (i in 1:length(Y_numcat)) {
      for (j in 1:nrow(Y)) {
        if (is.na(Y_cat[j,i])) {
          Yi[j,h:(h+Y_numcat[i]-2)]=NA
        }
      } 
      h=h+Y_numcat[i]-1
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
    betapost<- array(0, dim=c(nrow(betap),ncol(betap),(nimp-1)))
    bpost<-matrix(0,nrow(betap),ncol(betap))
    upost<-matrix(0,nrow(up),ncol(up))
    upostall<-array(0, dim=c(nrow(up),ncol(up),(nimp-1)))
    omegapost<- array(0, dim=c(nrow(covp),ncol(covp),(nimp-1)))
    opost<-matrix(0,nrow(covp),ncol(covp))
    covupost<- array(0, dim=c(nrow(covu),ncol(covu),(nimp-1)))
    cpost<-matrix(0,nrow(covu),ncol(covu))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    .Call("jomo1ranmix", Y, Yimp, Yimp2, Y_cat, X, Z, clus,betait,uit,bpost,upost,covit,opost, covuit, cpost, nburn, Sp,Sup,Y_numcat, 0,out.iter, PACKAGE = "jomo")
    #betapost[,,1]=bpost
    #upostall[,,1]=upost
    #omegapost[,,(1)]=opost
    #covupost[,,(1)]=cpost
    bpost<-matrix(0,nrow(betap),ncol(betap))
    opost<-matrix(0,nrow(covp),ncol(covp))
    upost<-matrix(0,nrow(up),ncol(up))
    cpost<-matrix(0,nrow(covu),ncol(covu))
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Y_cat
    if (output==1) cat("First imputation registered.", "\n")
    for (i in 2:nimp) {
      #Yimp2=matrix(0, nrow(Yimp),ncol(Yimp))
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
      imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+1):(ncol(Y)+ncol(X)+ncol(Z))]=Z
      imp[(i*nrow(clus)+1):((i+1)*nrow(clus)), (ncol(Y)+ncol(X)+ncol(Z)+1)]=clus
      imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+2)]=c(1:nrow(Y))
      imp[(i*nrow(Z)+1):((i+1)*nrow(Z)), (ncol(Y)+ncol(X)+ncol(Z)+3)]=i 
      .Call("jomo1ranmix", Y, Yimp, Yimp2, Y_cat, X, Z, clus,betait,uit,bpost,upost,covit,opost, covuit, cpost, nbetween, Sp,Sup,Y_numcat, 0, out.iter, PACKAGE = "jomo")
      betapost[,,(i-1)]=bpost
      upostall[,,(i-1)]=upost
      omegapost[,,(i-1)]=opost
      covupost[,,(i-1)]=cpost
      bpost<-matrix(0,nrow(betap),ncol(betap))
      opost<-matrix(0,nrow(covp),ncol(covp))
      upost<-matrix(0,nrow(up),ncol(up))
      cpost<-matrix(0,nrow(covu),ncol(covu))
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),1:ncol(Y)]=Y_cat
      if (output==1) cat("Imputation number ", i, "registered", "\n")
    }
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
    for (i in 1:ncol(Y)) {
      imp[,i]<-as.factor(imp[,i]) 
      levels(imp[,i])<-previous_levels[[i]]
    }
    if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y_cat), sep = "")
    if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycat,colnamx,colnamz,"clus","id","Imputation")
    return(imp)
  }
