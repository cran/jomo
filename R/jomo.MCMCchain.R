jomo.MCMCchain <-
  function(Y, X=NULL, Z=NULL,clus=NULL, betap=NULL, up=NULL, covp=NULL, covu=NULL, Sp=NULL, Sup=NULL, nburn=500, a=NULL, meth="common",output=1, out.iter=10) {
    if (is.null(clus)) {
      cat("No clustering, using functions for single level imputation.\n")
      imp<-jomo1.MCMCchain(Y, X, betap, covp, Sp, nburn, output, out.iter)
    }
    if (!is.null(clus)) {
      cat("Clustered data, using functions for two-level imputation.\n")
      imp<-jomo1ran.MCMCchain(Y, X, Z,clus, betap, up, covp, covu, Sp, Sup, nburn,a, meth, output, out.iter) 
    }
    return(imp)
  }
