\name{jomo.clmm.MCMCchain}
\alias{jomo.clmm.MCMCchain}

\title{
  clmm Compatible JM Imputation - A tool to check convergence of the MCMC
}
\description{
  This function is similar to the jomo.clmm function, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
  jomo.clmm.MCMCchain(formula, data, level=rep(1,ncol(data)), 
                      beta.start=NULL, l2.beta.start=NULL, u.start=NULL,
                      l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL,
                      l2cov.prior=NULL, a.start=NULL, a.prior=NULL,
                      betaY.start=NULL,  covuY.start=NULL,
                      uY.start=NULL, nburn=1000, meth="common", 
                      start.imp=NULL, start.imp.sub=NULL, l2.start.imp=NULL,
                      output=1, out.iter=10) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    an object of class formula: a symbolic description of the model to be fitted. It is possible to include in this formula interactions (through symbols '*' and '%') and polynomial terms (with the usual lm syntax, e.g. for a quadratic effect for variable x, 'I(x^2)'). Random effects follow the usual lmer syntax; however, it is possible only to allow for one grouping variable and all the variables with random effects must be included within the same brackets (es: (1+X|clus) is correct, while (1|clus)+(X|clus) is NOT).
  }
  
  \item{data}{
    A data.frame containing all the variables to include in the imputation model. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
  }
  \item{level}{
  A vector, indicating whether each variable is either a level 1 or a level 2 variable. The value assigned to the cluster indicator is irrelevant.   

  }
   \item{beta.start}{
    Starting value for beta, the vector(s) of level-1 fixed effects for the joint model for the covariates. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
  \item{l2.beta.start}{
    Starting value for beta2, the vector(s) of level-2 fixed effects for the joint model for the covariates. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
  \item{u.start}{
    A matrix where different rows are the starting values within each cluster of the random effects estimates u for the joint model for the covariates. The default is a matrix of zeros.
  }
  \item{l1cov.start}{
    Starting value of the level-1 covariance matrix of the joint model for the covariates. Dimension of this square matrix is equal to the number of covariates (continuous plus latent normals) in the imputation model. The default is the identity matrix. Functions for imputation with random cluster-specific covariance matrices are an exception, because we need to pass the starting values for all of the matrices stacked one above the other. 
  }
  \item{l2cov.start}{
    Starting value for the level 2 covariance matrix of the joint model for the covariates. Dimension of this square matrix is equal to the number of level-1 covariates (continuous plus latent normals) in the analysis model times the number of random effects plus the number of level-2 covariates. The default is an identity matrix.
  }
  \item{l1cov.prior}{
    Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
  }
  \item{l2cov.prior}{
    Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
  }
 \item{a.start}{
Starting value for the degrees of freedom of the inverse Wishart distribution of the cluster-specific covariance matrices. Default is 50+D, with D being the dimension of the covariance matrices. This is used only with clustered data and when option meth is set to "random".
}
  \item{a.prior}{
Hyperparameter (Degrees of freedom) of the chi square prior distribution for the degrees of freedom of the inverse Wishart distribution for the cluster-specific covariance matrices. Default is D, with D being the dimension of the covariance matrices.
}
  
  \item{meth}{
    Method used to deal with level 1 covariance matrix. When set to "common", a common matrix across clusters is used (functions jomo1rancon, jomo1rancat and jomo1ranmix). When set to "fixed", fixed study-specific matrices are considered (jomo1ranconhr, jomo1rancathr and jomo1ranmixhr with coption meth="fixed"). Finally, when set to "random", random study-specific matrices are considered (jomo1ranconhr, jomo1rancathr and jomo1ranmixhr with option meth="random")
  }

   \item{betaY.start}{
    Starting value for betaY, the vector of fixed effects for the substantive analysis model.  The default is the complete records estimate.
  }
  
   \item{covuY.start}{
    Starting value for covuY, the random effects covariance matrix of the substantive analysis model.  The default is the complete records estimate.
  }
  
    \item{uY.start}{
    Starting value for uY, the random effects matrix of the substantive analysis model.  The default is the complete records estimate.
  }
  
  \item{nburn}{
    Number of burn in iterations. Default is 1000.
  }
  \item{output}{
    When set to 0, no output is shown on screen at the end of the process. When set to 1, only the parameter estimates related to the substantive model are shown (default). When set to 2, all parameter estimates (posterior means) are displayed.
  }
  \item{out.iter}{
    When set to K, every K iterations a dot is printed on screen. Default is 10. 
  }
  \item{start.imp}{
    Starting value for the missing data in the covariates of the substantive model. n-level categorical variables are substituted by n-1 latent normals.
  }
    \item{l2.start.imp}{
    Starting value for the missing data in the level-2 covariates of the substantive model. n-level categorical variables are substituted by n-1 latent normals.
  }
  \item{start.imp.sub}{
    Starting value for the missing data in the outcome of the substantive model. For family="binomial", these are the values of the latent normals. 
  }

}

\value{
  A list is returned; this contains the final imputed dataset (finimp) and several 3-dimensional matrices, containing all the values drawn for each parameter at each iteration: these are fixed effect parameters of the covariates beta (collectbeta), level 1 covariance matrices (collectomega), fixed effect estimates of the substantive model and associated residual variances. If there are some categorical outcomes, a further output is included in the list, finimp.latnorm, containing the final state of the imputed dataset with the latent normal variables.   

}

\examples{

  # make sure social is a factor:
  
  cldata<-within(cldata, social<-factor(social))
  
  # we define the data frame with all the variables 
  
  data<-cldata[,c("measure","age", "social", "city")]
  
  # And the formula of the substantive lm model 
  # social as an outcome only because it is the only ordinal variable in the dataset...
  
  formula<-as.formula(social~age+measure+(1|city))
  
  #And finally we run the imputation function:
  
  imp<-jomo.clmm.MCMCchain(formula,data, nburn=100)
  
  # We can check, for example, the convergence of the first element of beta:
  
  # plot(c(1:100),imp$collectbeta[1,1,1:100],type="l")
  
  
}
