jomo1ran <-
  function(Y, X=NULL, Z=NULL,clus, beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, nburn=500, nbetween=100, nimp=5, a=NULL, meth="common",output=1, out.iter=10) {
    stopifnot(meth=="common"|meth=="fixed"|meth=="random")
    ncon=0
    ncat=0
    Y.con=NULL
    Y.cat=NULL
    Y.numcat=NULL
    for (i in 1:ncol(Y)) {
      if (is.numeric(Y[,i])) {
        ncon=ncon+1
        Y.con<-cbind(Y.con,Y[,i])
        colnames(Y.con)[ncon]<-colnames(Y)[i]
      }
      else {
        if (is.factor(Y[,i])) {
        ncat=ncat+1
        Y.cat<-cbind(Y.cat,Y[,i])
        colnames(Y.cat)[ncat]<-colnames(Y)[i]
        Y.numcat<-cbind(Y.numcat,nlevels(Y[,i]))
        }
      }
    }
    if (is.null(X)) X=matrix(1,nrow(Y),1)
    if (is.null(Z)) Z=matrix(1,nrow(Y),1)
    if (meth=="common") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1con.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),ncol(Y))
        if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(Y))
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(Y)*ncol(Z))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(Y))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(Y)*ncol(Z))
        imp<-jomo1rancon(Y.con, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp, output, out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancat.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),((sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(u.start))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(l2cov.start))
        imp<-jomo1rancat(Y.cat,Y.numcat, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp, output, out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmix.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(u.start))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(l2cov.start))
        imp<-jomo1ranmix(Y.con, Y.cat, Y.numcat, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp, output, out.iter)
      }
    }
    if (meth=="fixed") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1conhr with fixed cluster-l1cov.priorecific covariance matrices.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),ncol(Y))
        if (is.null(l1cov.start)) l1cov.start=matrix(diag(1,ncol(Y)),nrow(unique(clus))*ncol(Y),ncol(Y),2)
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(Y)*ncol(Z))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(Y))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(Y)*ncol(Z))
        imp<-jomo1ranconhr(Y.con, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp,a=15000,meth="fixed", output, out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with fixed cluster-l1cov.priorecific covariance matrices.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),((sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l1cov.start)) l1cov.start=matrix(diag(1,ncol(beta.start)),ncol(beta.start)*nrow(unique(clus)),ncol(beta.start),2)
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(u.start))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(l2cov.start))
        imp<-jomo1rancathr(Y.cat,Y.numcat, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp,a=15000, meth="fixed", output, out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmixhr with fixed cluster-l1cov.priorecific covariance matrices.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l1cov.start)) l1cov.start=l1cov.start=matrix(diag(1,ncol(beta.start)),ncol(beta.start)*nrow(unique(clus)),ncol(beta.start),2)
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(u.start))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(l2cov.start))
        imp<-jomo1ranmixhr(Y.con, Y.cat, Y.numcat, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp,a=15000, meth="fixed", output, out.iter)
      }
    }
    if (meth=="random") {
      if (ncat==0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomoran1conhr with random cluster-l1cov.priorecific covariance matrices.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),ncol(Y))
        if (is.null(l1cov.start)) l1cov.start=matrix(diag(1,ncol(Y)),nrow(unique(clus))*ncol(Y),ncol(Y),2)
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(Y)*ncol(Z))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(Y))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(Y)*ncol(Z))
        if (is.null(a)) a=ncol(l1cov.start)
        imp<-jomo1ranconhr(Y.con, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp,a=a,meth="random", output, out.iter)
      }
      if (ncat>0 & ncon==0) {
        cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with random cluster-l1cov.priorecific covariance matrices.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),((sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l1cov.start)) l1cov.start=matrix(diag(1,ncol(beta.start)),ncol(beta.start)*nrow(unique(clus)),ncol(beta.start),2)
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*((sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(u.start))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(l2cov.start))
        if (is.null(a)) a=ncol(l1cov.start)
        imp<-jomo1rancathr(Y.cat,Y.numcat, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp,a=a, meth="random", output, out.iter)
      }
      if (ncat>0 & ncon>0) {
        cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1ranmixhr with random cluster-l1cov.priorecific covariance matrices.", "\n")
        if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l1cov.start)) l1cov.start=l1cov.start=matrix(diag(1,ncol(beta.start)),ncol(beta.start)*nrow(unique(clus)),ncol(beta.start),2)
        if (is.null(u.start)) u.start=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
        if (is.null(l2cov.start)) l2cov.start=diag(1,ncol(u.start))
        if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(l1cov.start))
        if (is.null(l2cov.prior)) l2cov.prior=diag(1,ncol(l2cov.start))
        if (is.null(a)) a=ncol(l1cov.start)
        imp<-jomo1ranmixhr(Y.con, Y.cat, Y.numcat, X, Z, clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp,a=a, meth="random", output, out.iter)
      }
    }
    return(imp)
  }
