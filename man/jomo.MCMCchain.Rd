\name{jomo.MCMCchain}
\alias{jomo.MCMCchain}

\title{
  JM Imputation - A tool to check convergence of the MCMC
}
\description{
  This function is similar to the jomo function, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
  jomo.MCMCchain(Y, Y2=NULL, X=NULL, X2=NULL, Z=NULL, clus=NULL, 
                 beta.start=NULL, l2.beta.start=NULL, u.start=NULL,
                 l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, 
                l2cov.prior=NULL, start.imp=NULL, l2.start.imp=NULL, 
    nburn=1000, a=NULL, a.prior=NULL, meth="common",output=1, out.iter=10) 
  
}

\arguments{
  \item{Y}{
    A data.frame containing the outcomes of the imputation model, i.e. the partially observed level 1 variables. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
  }
   \item{Y2}{
    A data.frame containing the level-2 outcomes of the imputation model, i.e. the partially observed level-2 variables. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
  }
  \item{X}{
    A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
  }
  \item{X2}{
    A data frame, or matrix, with level-2 covariates of the joint imputation model. Rows correspond to different level-1 observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
  }
  \item{Z}{
    A data frame, or matrix, for covariates associated to random effects in the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
  }
  \item{clus}{
    A data frame, or matrix, containing the cluster indicator for each observation. If missing, functions for single level imputation are automatically used.
  }
  \item{beta.start}{
    Starting value for beta, the vector(s) of level-1 fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
  \item{l2.beta.start}{
    Starting value for beta2, the vector(s) of level-2 fixed effects. Rows index different covariates and columns index different level-2 outcomes. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
  \item{u.start}{
    A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
  }
  \item{l1cov.start}{
    Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix. Functions for imputation with random cluster-specific covariance matrices are an exception, because we need to pass the starting values for all of the matrices stacked one above the other. 
  }
  \item{l2cov.start}{
    Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model times the number of random effects plus the number of level-2 outcomes. The default is an identity matrix.
  }
  \item{l1cov.prior}{
    Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
  }
  \item{l2cov.prior}{
    Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
  }
   \item{start.imp}{
Starting value for the imputed dataset. n-level categorical variables are substituted by n-1 latent normals.
  }
   \item{l2.start.imp}{
Starting value for the level-2 imputed variables. n-level categorical variables are substituted by n-1 latent normals.
  }
  \item{nburn}{
    Number of iterations. Default is 1000.
  }
  \item{a}{
Starting value for the degrees of freedom of the inverse Wishart distribution of the cluster-specific covariance matrices. Default is 50+D, with D being the dimension of the covariance matrices. This is used only with clustered data and when option meth is set to "random".
}
  \item{a.prior}{
Hyperparameter (Degrees of freedom) of the chi square prior distribution for the degrees of freedom of the inverse Wishart distribution for the cluster-specific covariance matrices. Default is D, with D being the dimension of the covariance matrices.
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
  A list is returned; this contains the final imputed dataset (finimp) and several 3-dimensional matrices, containing all the values drawn for each parameter at each iteration: these are, potentially, fixed effect parameters beta (collectbeta), random effects (collectu), level 1 (collectomega) and level 2 covariance matrices (collectcovu) and level-2 fixed effect parameters. If there are some categorical outcomes, a further output is included in the list, finimp.latnorm, containing the final state of the imputed dataset with the latent normal variables.
}

\examples{
  #First of all we load and attach the data:
  
  data(cldata)
  attach(cldata)
  
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
