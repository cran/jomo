jomo <-
  function(Y, X=NULL, Z=NULL,clus=NULL, beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL, nburn=500, nbetween=100, nimp=5, a=NULL, meth="common",output=1, out.iter=10) {
    if (is.null(clus)) {
      cat("No clustering, using functions for single level imputation.\n")
      imp<-jomo1(Y, X, beta.start, l1cov.start, l1cov.prior, nburn, nbetween, nimp, output, out.iter)
    }
    if (!is.null(clus)) {
      cat("Clustered data, using functions for two-level imputation.\n")
      imp<-jomo1ran(Y, X, Z,clus, beta.start, u.start, l1cov.start, l2cov.start, l1cov.prior, l2cov.prior, nburn, nbetween, nimp, a, meth, output, out.iter) 
    }
    return(imp)
  }
