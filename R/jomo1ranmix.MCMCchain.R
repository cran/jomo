jomo1ranmix.MCMCchain <-
  function(Y_con, Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), Z=matrix(1,nrow(Y_cat),1), clus, betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))), up=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))), covp=diag(1,ncol(betap)), covu=diag(1,ncol(up)), Sp=diag(1,ncol(covp)), Sup=diag(1,ncol(covu)), nburn=100, output=1, out.iter=10) {
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
    stopifnot(nrow(Y_con)==nrow(clus),nrow(Y_con)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))),nrow(covp)==ncol(covp), nrow(covp)==ncol(betap), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp),nrow(Z)==nrow(Y_con), ncol(covu)==ncol(up), ncol(up)==ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
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
    nimp=1
    colnamycon<-colnames(Y_con)
    colnamycat<-colnames(Y_cat)
    colnamx<-colnames(X)
    colnamz<-colnames(Z)
    Y_con<-as.matrix(Y_con,nrow(Y_con),ncol(Y_con))
    Y_cat<-as.matrix(Y_cat,nrow(Y_cat),ncol(Y_cat))
    X<-as.matrix(X,nrow(X),ncol(X))
    Z<-as.matrix(Z,nrow(Z),ncol(Z))
    clus<-as.matrix(clus,nrow(clus),ncol(clus))
    Y=cbind(Y_con,Y_cat)
    Yi=cbind(Y_con, matrix(0,nrow(Y_con),(sum(Y_numcat)-length(Y_numcat))))
    h=1
    for (i in 1:length(Y_numcat)) {
      for (j in 1:nrow(Y)) {
        if (is.na(Y_cat[j,i])) {
          Yi[j,(ncol(Y_con)+h):(ncol(Y_con)+h+Y_numcat[i]-2)]=NA
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
    .Call("MCMCjomo1ranmix", Y, Yimp, Yimp2, Y_cat, X, Z, clus,betait,uit,betapost,upostall,covit,omegapost, covuit, covupost, nburn, Sp,Sup,Y_numcat, ncol(Y_con),out.iter, PACKAGE = "jomo")
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y_con)]=Yimp2[,1:ncol(Y_con)]
    imp[(nrow(Y)+1):(2*nrow(Y)),(ncol(Y_con)+1):ncol(Y)]=Y_cat
    for (i in 1:ncol(Y_cat)) {
      if (Ycatsum1[i]==1) {
        imp[,(ncol(Y_con)+i)]<-factor(as.numeric(imp[,(ncol(Y_con)+i)])-1)                   
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
    if (is.null(colnamycon)) colnamycon=paste("Ycon", 1:ncol(Y_con), sep = "")
    if (is.null(colnamz)) colnamz=paste("Z", 1:ncol(Z), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycon,colnamycat,colnamx,colnamz,"clus","id","Imputation")
    return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost,"collectu"=upostall, "collectcovu"=covupost))
  }
