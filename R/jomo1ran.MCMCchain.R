jomo1ran.MCMCchain <-
  function(Y, X=NULL, Z=NULL,clus, betap=NULL, up=NULL, covp=NULL, covu=NULL, Sp=NULL, Sup=NULL, nburn=500, a=NULL, meth="common", output=1, out.iter=10) {
    stopifnot(meth=="common"|meth=="fixed"|meth=="random")
    ncon=0
    ncat=0
    Y_con=NULL
    Y_cat=NULL
    Y_numcat=NULL
    for (i in 1:ncol(Y)) {
      if (is.numeric(Y[,i])) {
        ncon=ncon+1
        Y_con<-cbind(Y_con,Y[,i])
        colnames(Y_con)[ncon]<-colnames(Y)[i]
      }
      else if (is.factor(Y[,i])) {
        ncat=ncat+1
        Y_cat<-cbind(Y_cat,Y[,i])
        colnames(Y_cat)[ncat]<-colnames(Y)[i]
        Y_numcat<-cbind(Y_numcat,max(as.numeric(Y[!is.na(Y[,2]),i])))
      }
    }
    if (is.null(X)) X=matrix(1,nrow(Y),1)
    if (is.null(Z)) Z=matrix(1,nrow(Y),1)
    if (meth=="common") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1con.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),ncol(Y))
        if (is.null(covp)) covp=diag(1,ncol(Y))
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y))
        if (is.null(covu)) covu=diag(1,ncol(Y)*ncol(Z))
        if (is.null(Sp)) Sp=diag(1,ncol(Y))
        if (is.null(Sup)) Sup=diag(1,ncol(Y)*ncol(Z))
        imp<-jomo1rancon.MCMCchain(Y_con, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, output,out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancat.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covp)) covp=diag(1,ncol(betap))
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covu)) covu=diag(1,ncol(up))
        if (is.null(Sp)) Sp=diag(1,ncol(covp))
        if (is.null(Sup)) Sup=diag(1,ncol(covu))
        imp<-jomo1rancat.MCMCchain(Y_cat,Y_numcat, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, output, out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmix.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covp)) covp=diag(1,ncol(betap))
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covu)) covu=diag(1,ncol(up))
        if (is.null(Sp)) Sp=diag(1,ncol(covp))
        if (is.null(Sup)) Sup=diag(1,ncol(covu))
        imp<-jomo1ranmix.MCMCchain(Y_con, Y_cat, Y_numcat, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, output, out.iter)
      }
    }
    if (meth=="fixed") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1conhr with fixed cluster-specific covariance matrices.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),ncol(Y))
        if (is.null(covp)) covp=matrix(diag(1,ncol(Y)),nrow(unique(clus))*ncol(Y),ncol(Y),2)
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y))
        if (is.null(covu)) covu=diag(1,ncol(Y)*ncol(Z))
        if (is.null(Sp)) Sp=diag(1,ncol(Y))
        if (is.null(Sup)) Sup=diag(1,ncol(Y)*ncol(Z))
        imp<-jomo1ranconhr.MCMCchain(Y_con, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn,a=15000,meth="fixed", output, out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with fixed cluster-specific covariance matrices.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covp)) covp=matrix(diag(1,ncol(betap)),ncol(betap)*nrow(unique(clus)),ncol(betap),2)
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covu)) covu=diag(1,ncol(up))
        if (is.null(Sp)) Sp=diag(1,ncol(covp))
        if (is.null(Sup)) Sup=diag(1,ncol(covu))
        imp<-jomo1rancathr.MCMCchain(Y_cat,Y_numcat, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, a=15000, meth="fixed", output, out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmixhr with fixed cluster-specific covariance matrices.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covp)) covp=covp=matrix(diag(1,ncol(betap)),ncol(betap)*nrow(unique(clus)),ncol(betap),2)
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covu)) covu=diag(1,ncol(up))
        if (is.null(Sp)) Sp=diag(1,ncol(covp))
        if (is.null(Sup)) Sup=diag(1,ncol(covu))
        imp<-jomo1ranmixhr.MCMCchain(Y_con, Y_cat, Y_numcat, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, a=15000, meth="fixed", output, out.iter)
      }
    }
    if (meth=="random") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1conhr with random cluster-specific covariance matrices.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),ncol(Y))
        if (is.null(covp)) covp=matrix(diag(1,ncol(Y)),nrow(unique(clus))*ncol(Y),ncol(Y),2)
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y))
        if (is.null(covu)) covu=diag(1,ncol(Y)*ncol(Z))
        if (is.null(Sp)) Sp=diag(1,ncol(Y))
        if (is.null(Sup)) Sup=diag(1,ncol(Y)*ncol(Z))
        if (is.null(a)) a=ncol(covp)
        imp<-jomo1ranconhr.MCMCchain(Y_con, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, a,meth="random", output, out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with random cluster-specific covariance matrices.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covp)) covp=matrix(diag(1,ncol(betap)),ncol(betap)*nrow(unique(clus)),ncol(betap),2)
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covu)) covu=diag(1,ncol(up))
        if (is.null(Sp)) Sp=diag(1,ncol(covp))
        if (is.null(Sup)) Sup=diag(1,ncol(covu))
        if (is.null(a)) a=ncol(covp)
        imp<-jomo1rancathr.MCMCchain(Y_cat,Y_numcat, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, a, meth="random", output, out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmixhr with random cluster-specific covariance matrices.", "\n")
        if (is.null(betap)) betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covp)) covp=covp=matrix(diag(1,ncol(betap)),ncol(betap)*nrow(unique(clus)),ncol(betap),2)
        if (is.null(up)) up=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
        if (is.null(covu)) covu=diag(1,ncol(up))
        if (is.null(Sp)) Sp=diag(1,ncol(covp))
        if (is.null(Sup)) Sup=diag(1,ncol(covu))
        if (is.null(a)) a=ncol(covp)
        imp<-jomo1ranmixhr.MCMCchain(Y_con, Y_cat, Y_numcat, X, Z, clus, betap, up, covp, covu, Sp, Sup, nburn, a, meth="random", output, out.iter)
      }
    }
    return(imp)
  }
