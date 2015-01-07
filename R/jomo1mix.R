jomo1mix <-
  function(Y_con, Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))), covp=diag(1,ncol(betap)), Sp=diag(1,ncol(betap)), nburn=100, nbetween=100, nimp=5, meth="MH") {
    stopifnot((meth=="MH"|meth=="IW"),nrow(Y_con)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))),nrow(covp)==ncol(covp), nrow(covp)==ncol(betap), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp))
    rngflag=0
    colnamycon<-colnames(Y_con)
    colnamycat<-colnames(Y_cat)
    colnamx<-colnames(X)
    Y_con<-as.matrix(Y_con,nrow(Y_con),ncol(Y_con))
    Y_cat<-as.matrix(Y_cat,nrow(Y_cat),ncol(Y_cat))
    X<-as.matrix(X,nrow(X),ncol(X))
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
    if (meth=="MH") {
      .Call("jomo1mix", Y, Yimp, Yimp2, Y_cat, X,betap,bpost,covp,opost, nburn, Sp,Y_numcat, ncol(Y_con),rngflag, PACKAGE = "jomo")
    }
    if (meth=="IW") {
      .Call("jomo1mix2", Y, Yimp, Yimp2, Y_cat, X,betap,bpost,covp,opost, nburn, Sp,Y_numcat, ncol(Y_con),rngflag, PACKAGE = "jomo")
    }
    #betapost[,,1]=bpost
    #omegapost[,,1]=opost
    bpost<-matrix(0,nrow(betap),ncol(betap))
    opost<-matrix(0,nrow(covp),ncol(covp))
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y_con)]=Yimp2[,1:ncol(Y_con)]
    imp[(nrow(Y)+1):(2*nrow(Y)),(ncol(Y_con)+1):ncol(Y)]=Y_cat
    cat("First imputation registered.", "\n")
    for (i in 2:nimp) {
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
      imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+1)]=i
      imp[(i*nrow(X)+1):((i+1)*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
      if (meth=="MH") {
        .Call("jomo1mix", Y, Yimp, Yimp2, Y_cat, X,betap,bpost,covp, opost, nbetween, Sp, Y_numcat, ncol(Y_con),rngflag, PACKAGE = "jomo") 
      }
      if (meth=="IW") {
        .Call("jomo1mix2", Y, Yimp, Yimp2, Y_cat, X,betap,bpost,covp, opost, nbetween, Sp, Y_numcat, ncol(Y_con),rngflag, PACKAGE = "jomo") 
      }
      betapost[,,(i-1)]=bpost
      omegapost[,,(i-1)]=opost
      bpost<-matrix(0,nrow(betap),ncol(betap))
      opost<-matrix(0,nrow(covp),ncol(covp))
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),1:ncol(Y_con)]=Yimp2[,1:ncol(Y_con)]
      imp[(i*nrow(X)+1):((i+1)*nrow(X)),(ncol(Y_con)+1):ncol(Y)]=Y_cat
      cat("Imputation number ", i, "registered", "\n")
    }
    betapostmean<-apply(betapost, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    cat("The posterior mean of the fixed effects estimates is:\n")
    print(betapostmean)
    cat("The posterior covariance matrix is:\n")
    print(omegapostmean)
    imp<-data.frame(imp)
    if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y_cat), sep = "")
    if (is.null(colnamycon)) colnamycon=paste("Ycon", 1:ncol(Y_con), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycon,colnamycat,colnamx,"Imputation","id")
    return(imp)
  }
