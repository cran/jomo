jomo1mix.MCMCchain <-
  function(Y_con, Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))), covp=diag(1,ncol(betap)), Sp=diag(1,ncol(betap)), nburn=100, output=1, out.iter=10) {
    previous_levels<-list()
    for (i in 1:ncol(Y_cat)) {
      Y_cat[,i]<-factor(Y_cat[,i])
      previous_levels[[i]]<-levels(Y_cat[,i])
      levels(Y_cat[,i])<-1:nlevels(Y_cat[,i])
    }
    for (i in 1:ncol(X)) {
      if (is.factor(X[,i])) X[,i]<-as.numeric(X[,i])
    }
    stopifnot(nrow(Y_con)==nrow(X), nrow(betap)==ncol(X), ncol(betap)==(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))),nrow(covp)==ncol(covp), nrow(covp)==ncol(betap), nrow(Sp)==ncol(Sp),nrow(Sp)==nrow(covp))
    betait=matrix(0,nrow(betap),ncol(betap))
    for (i in 1:nrow(betap)) {
      for (j in 1:ncol(betap)) betait[i,j]=betap[i,j]
    }
    covit=matrix(0,nrow(covp),ncol(covp))
    for (i in 1:nrow(covp)) {
      for (j in 1:ncol(covp)) covit[i,j]=covp[i,j]
    }   
    nimp=1
    colnamycon<-colnames(Y_con)
    colnamycat<-colnames(Y_cat)
    colnamx<-colnames(X)
    Y_con<-data.matrix(Y_con)
    Y_cat<-data.matrix(Y_cat)
    X<-data.matrix(X)
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
    imp=matrix(0,nrow(Y)*(nimp+1),ncol(Y)+ncol(X)+2)
    imp[1:nrow(Y),1:ncol(Y)]=Y
    imp[1:nrow(X), (ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[1:nrow(X), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
    Yimp=Yi
    Yimp2=matrix(Yimp, nrow(Yimp),ncol(Yimp))
    imp[(nrow(X)+1):(2*nrow(X)),(ncol(Y)+1):(ncol(Y)+ncol(X))]=X
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+1)]=1
    imp[(nrow(X)+1):(2*nrow(X)), (ncol(Y)+ncol(X)+2)]=c(1:nrow(Y))
    betapost<- array(0, dim=c(nrow(betap),ncol(betap),nburn))
    omegapost<- array(0, dim=c(nrow(covp),ncol(covp),nburn))
    meanobs<-colMeans(Yi,na.rm=TRUE)
    for (i in 1:nrow(Yi)) for (j in 1:ncol(Yi)) if (is.na(Yimp[i,j])) Yimp2[i,j]=meanobs[j]
    .Call("MCMCjomo1mix", Y, Yimp, Yimp2, Y_cat, X,betait,betapost,covit,omegapost, nburn, Sp,Y_numcat, ncol(Y_con),out.iter, PACKAGE = "jomo")
    imp[(nrow(Y)+1):(2*nrow(Y)),1:ncol(Y_con)]=Yimp2[,1:ncol(Y_con)]
    imp[(nrow(Y)+1):(2*nrow(Y)),(ncol(Y_con)+1):ncol(Y)]=Y_cat
    betapostmean<-apply(betapost, c(1,2), mean)
    omegapostmean<-apply(omegapost, c(1,2), mean)
    if (output==1) {
      cat("The posterior mean of the fixed effects estimates is:\n")
      print(betapostmean)
      cat("The posterior covariance matrix is:\n")
      print(omegapostmean)
    }
    imp<-data.frame(imp)
    for (i in 1:ncol(Y_cat)) {
      imp[,(ncol(Y_con)+i)]<-as.factor(imp[,(ncol(Y_con)+i)]) 
      levels(imp[,(ncol(Y_con)+i)])<-previous_levels[[i]]
    }
    if (is.null(colnamycat)) colnamycat=paste("Ycat", 1:ncol(Y_cat), sep = "")
    if (is.null(colnamycon)) colnamycon=paste("Ycon", 1:ncol(Y_con), sep = "")
    if (is.null(colnamx)) colnamx=paste("X", 1:ncol(X), sep = "")
    colnames(imp)<-c(colnamycon,colnamycat,colnamx,"Imputation","id")
    return(list("finimp"=imp,"collectbeta"=betapost,"collectomega"=omegapost))
  }
