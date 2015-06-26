jomo1rancathr.MCMCchain <-
  function(Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), Z=matrix(1,nrow(Y_cat),1), clus, betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat)))), up=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y_numcat)-length(Y_numcat)))), covp=matrix(diag(1,ncol(betap)),ncol(betap)*nrow(unique(clus)),ncol(betap),2), covu=diag(1,ncol(up)), Sp=diag(1,ncol(betap)), Sup=diag(1,ncol(up)),nburn=100,a=ncol(betap),meth="random", output=1, out.iter=10) {
    Ycatsum1<-rep(0,ncol(Y_cat))
    for (i in 1:ncol(Y_cat)) {
      if (min(as.numeric(Y_cat[!is.na(Y_cat[,i]),i]))==0) {
        Y_cat[,i]<-factor(as.numeric(Y_cat[,i])+1)
        Ycatsum1[i]<-1
      }
    }
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    for (i in 1:ncol(Z)) {
      if (is.factor(Z[,i])) Z[,i]<-as.numeric(Z[,i])
    }
    stopifnot((meth=="fixed"|meth=="random"),nrow(betap)==ncol(X), ncol(betap)==((sum(Y_numcat)-length(Y_numcat))),nrow(covp)==nrow(up)*ncol(covp), nrow(covp)==nrow(up)*ncol(betap), nrow(Sp)==ncol(Sp),nrow(covp)==nrow(up)*nrow(Sp),nrow(Z)==nrow(Y_cat), ncol(covu)==ncol(up), ncol(up)==ncol(Z)*((sum(Y_numcat)-length(Y_numcat))),det(Sp)>0)
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
    colnamycat<-colnames(Y_cat)
    colnamx<-colnames(X)
    colnamz<-colnames(Z)
    Y_cat<-as.matrix(Y_cat,nrow(Y_cat),ncol(Y_cat))
    X<-as.matrix(X,nrow(X),ncol(X))
    Z<-as.matrix(Z,nrow(Z),ncol(Z))
    clus<-as.matrix(clus,nrow(clus),ncol(clus))
    Y=cbind(Y_cat)
    Yi=cbind( matrix(0,nrow(Y_cat),(sum(Y_numcat)-length(Y_numcat))))
    h=1
    for (i in 1:length(Y_numcat)) {
      for (j in 1:nrow(Y)) {
        if (is.na(Y_cat[j,i])) {
          Yi[j,(h):(h+Y_numcat[i]-2)]=NA
        }
      } 
      h=h+Y_numcat[i]-1
    }
    if (output!=1) out.iter=nburn+2
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
    betapost<- array(0, dim=c(nrow(betap),ncol(betap),nburn))
    omegapost<- array(0, dim=c(nrow(covp),ncol(covp),nburn))
    upostall<-array(0, dim=c(nrow(up),ncol(up),nburn))
    covupost<- array(0, dim=c(nrow(covu),ncol(covu),nburn))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    if (meth=="fixed") {
      .Call("MCMCjomo1ranmixhf", Y, Yimp, Yimp2, Y_cat, X, Z, clus,betait,uit,betapost,upostall,covit,omegapost, covuit,covupost,nburn, Sp,Sup,Y_numcat, 0,ait,out.iter, PACKAGE = "jomo")
    }
    if (meth=="random") {
      .Call("MCMCjomo1ranmixhr", Y, Yimp, Yimp2, Y_cat, X, Z, clus,betait,uit,betapost,upostall,covit,omegapost, covuit,covupost,nburn, Sp,Sup,Y_numcat, 0,ait,out.iter, PACKAGE = "jomo")
    }
    
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Y_cat
    for (i in 1:ncol(Y)) {
      if (Ycatsum1[i]==1) {
        imp[,i]<-factor(as.numeric(imp[,i])-1)                   
      }
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
    if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y_cat), sep = "")
    if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycat,colnamx,colnamz,"clus","id","Imputation")
    return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost,"collectu"=upostall, "collectcovu"=covupost))
  }
