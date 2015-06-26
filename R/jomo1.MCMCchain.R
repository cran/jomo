jomo1.MCMCchain <-
  function(Y, X=NULL, betap=NULL, covp=NULL, Sp=NULL, nburn=500, output=1, out.iter=10) {
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
        Y_numcat<-cbind(Y_numcat,max(as.numeric(Y[!is.na(Y[,i]),i])))
      }
    }
    if (is.null(X)) X=matrix(1,nrow(Y),1)
    if (ncat==0 & ncon>0) {
      cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomo1con.", "\n")
      if (is.null(betap)) betap=matrix(0,ncol(X),(ncol(Y_con)))
      if (is.null(covp)) covp=diag(1,ncol(betap))
      if (is.null(Sp)) Sp=diag(1,ncol(betap))
      imp<-jomo1con.MCMCchain(Y_con, X, betap, covp, Sp, nburn, output,out.iter)
    }
    if (ncat>0 & ncon==0) {
      cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1cat.", "\n")
      if (is.null(betap)) betap=matrix(0,ncol(X),(sum(Y_numcat)-length(Y_numcat)))
      if (is.null(covp)) covp=diag(1,ncol(betap))
      if (is.null(Sp)) Sp=diag(1,ncol(betap))
      imp<-jomo1cat.MCMCchain(Y_cat,Y_numcat,X, betap, covp, Sp, nburn, output,out.iter)
    }
    if (ncat>0 & ncon>0) {
      cat("Found ", ncon, "continuous outcomes and ", ncat, "categorical. Using function jomo1mix.", "\n")
      if (is.null(betap)) betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat))))
      if (is.null(covp)) covp=diag(1,ncol(betap))
      if (is.null(Sp)) Sp=diag(1,ncol(betap))
      imp<-jomo1mix.MCMCchain(Y_con,Y_cat,Y_numcat,X, betap, covp, Sp, nburn, output,out.iter)
    }
    return(imp)
  }
