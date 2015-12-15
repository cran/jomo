jomo1.MCMCchain <-
  function(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, nburn=500, output=1, out.iter=10) {
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
    if (ncat==0 & ncon>0) {
      cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomo1con.", "\n")
      if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(ncol(Y.con)))
      if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
      if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
      imp<-jomo1con.MCMCchain(Y.con, X, beta.start, l1cov.start, l1cov.prior, nburn, output,out.iter)
    }
    if (ncat>0 & ncon==0) {
      cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1cat.", "\n")
      if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(sum(Y.numcat)-length(Y.numcat)))
      if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
      if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
      imp<-jomo1cat.MCMCchain(Y.cat,Y.numcat,X, beta.start, l1cov.start, l1cov.prior, nburn, output,out.iter)
    }
    if (ncat>0 & ncon>0) {
      cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1mix.", "\n")
      if (is.null(beta.start)) beta.start=matrix(0,ncol(X),(ncol(Y.con)+(sum(Y.numcat)-length(Y.numcat))))
      if (is.null(l1cov.start)) l1cov.start=diag(1,ncol(beta.start))
      if (is.null(l1cov.prior)) l1cov.prior=diag(1,ncol(beta.start))
      imp<-jomo1mix.MCMCchain(Y.con,Y.cat,Y.numcat,X, beta.start, l1cov.start, l1cov.prior, nburn, output,out.iter)
    }
    return(imp)
  }
