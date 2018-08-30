\name{jomo.coxph.MCMCchain}
\alias{jomo.coxph.MCMCchain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 coxph Compatible JM Imputation - A tool to check convergence of the MCMC
}
\description{
 This function is similar to the jomo.coxph function, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
  jomo.coxph.MCMCchain(formula, data, beta.start = NULL, l1cov.start = NULL,
                 l1cov.prior = NULL, nburn = 1000, start.imp = NULL,
                 betaY.start = NULL, output = 1, out.iter = 10)
  }
\arguments{
    \item{formula}{
      an object of class formula: a symbolic description of the model to be fitted. It is possible to include in this formula interactions (through symbols '*' and '%') and polynomial terms (with the usual lm syntax, e.g. for a quadratic effect for variable x, 'I(x^2)'). Survival model is specified with the usual coxph syntax.
    }
    \item{data}{
      A data.frame containing all the variables to include in the imputation model. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
    }
    
    \item{beta.start}{
      Starting value for beta, the vector(s) of fixed effects for the joint model for the covariates. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
    }
    
    \item{l1cov.start}{
      Starting value of the level-1 covariance matrix of the joint model for the covariates. Dimension of this square matrix is equal to the number of covariates (continuous plus latent normals) in the imputation model. The default is the identity matrix. 
    }
    
    \item{l1cov.prior}{
      Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
    }
     \item{betaY.start}{
    Starting value for betaY, the vector of fixed effects for the substantive analysis model.  The default is the complete records estimate.
  }
    \item{nburn}{
      Number of burn in iterations. Default is 1000.
    }
    \item{output}{
      When set to any value different from 1 (default), no output is shown on screen at the end of the process.
    }
    \item{out.iter}{
      When set to K, every K iterations a dot is printed on screen. Default is 10. 
    }
      \item{start.imp}{
Starting value for the missing data in the covariates of the substantive model. n-level categorical variables are substituted by n-1 latent normals.
  }
    
  }

  
  \value{
 A list is returned; this contains the final imputed dataset (finimp) and several 3-dimensional matrices, containing all the values drawn for each parameter at each iteration: these are fixed effect parameters of the covariates beta (collectbeta), level 1 covariance matrices (collectomega), fixed effect estimates of the substantive model. If there are some categorical outcomes, a further output is included in the list, finimp.latnorm, containing the final state of the imputed dataset with the latent normal variables.
  }

  \examples{
    # define substantive model

    formula<-as.formula(Surv(time, status) ~ measure + sex + I(measure^2))
    
    #Run imputation
    
if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  imp<-jomo.coxph.MCMCchain(formula,surdata, nburn = 100)
  }
    
}
  
  