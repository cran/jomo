\name{jomo.MCMCchain}
\alias{jomo.MCMCchain}

\title{
  JM Imputation - A tool to check convergence of the MCMC
}
\description{
  This function is similar to the jomo function, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
  jomo.MCMCchain(Y, X=NULL, Z=NULL,clus=NULL, beta.start=NULL, u.start=NULL,
                     l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, 
                     l2cov.prior=NULL, nburn=500, a=NULL, meth="common", 
                     output=1, out.iter=10) 
  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
    a data.frame where columns related to continuous variables are numeric and columns related to binary/categorical variables are factors. 
  }
  \item{X}{
    A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
  }
  \item{Z}{
    A data frame, or matrix, for covariates associated to random effects in the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
  }
  \item{clus}{
    A data frame, or matrix, containing the cluster indicator for each observation. If missing, functions for single level imputation are automatically used.
  }
  \item{beta.start}{
    Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we define n-1 latent normals. The default is a matrix of zeros.
  }
  \item{u.start}{
    A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
  }
  \item{l1cov.start}{
    Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
  }
  \item{l2cov.start}{
    Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model times the number of random effects. The default is an identity matrix.
  }
  \item{l1cov.prior}{
    Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
  }
  \item{l2cov.prior}{
    Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
  }
  \item{nburn}{
    Number of iterations. Default is 100.
  }
  \item{a}{
    Starting value for the degrees of freedom of the inverse Wishart distribution from which all of the covariance matrices are drawn in case method="random" option is specified. Default is the minimum possible, i.e. the dimension of the covariance matrices.
  }
  \item{meth}{
    Method used to deal with level 1 covariance matrix. When set to "common", a common matrix across clusters is used (functions jomo1rancon, jomo1rancat and jomo1ranmix). When set to "fixed", fixed study-specific matrices are considered (jomo1ranconhr, jomo1rancathr and jomo1ranmixhr with coption meth="fixed"). Finally, when set to "random", random study-specific matrices are considered (jomo1ranconhr, jomo1rancathr and jomo1ranmixhr with option meth="random")
  }
  
  \item{output}{
    When set to any value different from 1 (default), no output is shown on screen at the end of the process.
  }
  \item{out.iter}{
    When set to K, every K iterations a message "Iteration number N*K completed" is printed on screen. Default is 10. 
  }
  
}

\value{
  A list with five elements is returned: the final imputed dataset (finimp) and four 3-dimensional matrices, containing all the values for beta (collectbeta), the random effects (collectu) and the level 1 (collectomega) and level 2 covariance matrices (collectcovu). If functions for single level imputation are used, no cllectu or collectcovu are present.
}

\examples{
  #First of all we load and atach the data:
  
  data(mldata)
  attach(mldata)
  
  #Then we define all the inputs:
  
  Y<-data.frame(measure,age)
  clus<-data.frame(city)
  nburn=as.integer(200);
  
  #And finally we run the imputation function:
  imp<-jomo.MCMCchain(Y,clus=clus,nburn=nburn)
  #We can check the convergence of the first element of beta:
  
  plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")
  
  #Or similarly we can check the convergence of any element of th elevel 2 covariance matrix:
  
  plot(c(1:nburn),imp$collectcovu[1,2,1:nburn],type="l")
}
