jomo1cat <-
  function(Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat)))), covp=diag(1,ncol(betap)), Sp=diag(1,ncol(betap)), nburn=100, nbetween=100, nimp=5, output=1, out.iter=10) {
    previous_levels<-list()
    for (i in 1:ncol(Y_cat)) {
      Y_cat[,i]<-factor(Y_cat[,i])
      previous_levels[[i]]<-levels(Y_cat[,i])
      levels(Y_cat[,i])<-1:nlevels(Y_cat[,i])
    }
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    stopifnot( nrow(betap)==ncol(X), ncol(betap)==((sum(Y_numcat)-length(Y_numcat))),nrow(covp)==ncol(covp), nrow(covp)==ncol(betap), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp))
    betait=matrix(0,nrow(betap),ncol(betap))
    for (i in 1:nrow(betap)) {
      for (j in 1:ncol(betap)) betait[i,j]=betap[i,j]
    }
    covit=matrix(0,nrow(covp),ncol(covp))
    for (i in 1:nrow(covp)) {
      for (j in 1:ncol(covp)) covit[i,j]=covp[i,j]
    }    
    colnamycat<-colnames(Y_cat)
    colnamx<-colnames(X)
    Y_cat<-data.matrix(Y_cat)
    storage.mode(Y_cat) <- "numeric"    
    X<-data.matrix(X)
    storage.mode(X) <- "numeric"
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
    imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+2)
    imp[1:nrow(Y),1:ncol(Y)]=Y
    imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[1:nrow(X), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
    Yimp=Yi
    Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
    imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+1)]=1
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
    betapost<- array(0, dim=c(nrow(betap),ncol(betap),(nimp-1)))
    bpost<-matrix(0,nrow(betap),ncol(betap))
    omegapost<- array(0, dim=c(nrow(covp),ncol(covp),(nimp-1)))
    opost<-matrix(0,nrow(covp),ncol(covp))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    .Call("jomo1mix", Y, Yimp, Yimp2, Y_cat, X,betait,bpost,covit,opost, nburn, Sp,Y_numcat, 0, out.iter, PACKAGE = "jomo")
    #betapost[,,1]=bpost
    #omegapost[,,1]=opost
    bpost<-matrix(0,nrow(betap),ncol(betap))
    opost<-matrix(0,nrow(covp),ncol(covp))
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y)]=Y_cat
    if (output==1) cat("First imputation registered.", "\n")
    for (i in 2:nimp) {
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
      imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+1)]=i
      imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
      .Call("jomo1mix", Y, Yimp, Yimp2, Y_cat, X,betait,bpost,covit, opost, nbetween, Sp, Y_numcat, 0, out.iter, PACKAGE = "jomo") 
      
      betapost[,,(i-1)]=bpost
      omegapost[,,(i-1)]=opost
      bpost<-matrix(0,nrow(betap),ncol(betap))
      opost<-matrix(0,nrow(covp),ncol(covp))      
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),1:ncol(Y)]=Y_cat
      if (output==1) cat("Imputation number ", i, "registered", "\n")
    }

    betapostmean<-apply(betapost, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    if (output==1) {
      cat("The posterior mean of the fixed effects estimates is:\n")
      print(betapostmean)
      cat("The posterior covariance matrix is:\n")
      print(omegapostmean)
    }
    imp<-data.frame(imp)
    for (i in 1:ncol(Y)) {
        imp[,i]<-as.factor(imp[,i]) 
        levels(imp[,i])<-previous_levels[[i]]
    }
    if (is.null(colnamycat)) colnamycat=paste("Y", 1:ncol(Y_cat), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycat,colnamx,"Imputation","id")
    return(imp)
  }
