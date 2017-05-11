\name{jomo.glm.MCMCchain}
\alias{jomo.glm.MCMCchain}

\title{
  glm Compatible JM Imputation - A tool to check convergence of the MCMC
}
\description{
  This function is similar to the jomo.glm function, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
  jomo.glm.MCMCchain(formula, data, start.imp=NULL, start.imp.sub=NULL, 
                    beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL,
                    betaY.start=NULL, varY.start=NULL, nburn=1000, 
                    output=1, out.iter=10, family="binomial") 
  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{formula}{
    an object of class formula: a symbolic description of the model to be fitted. It is possible to include in this formula interactions (through symbols '*' and '%') and polynomial terms (with the usual glm syntax, e.g. for a quadratic effect for variable x, 'I(x^2)'). 
  }
  
   \item{data}{
    A data.frame containing all the variables to include in the imputation model. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
  }
   \item{start.imp}{
Starting value for the imputed covariates. n-level categorical variables are substituted by n-1 latent normals.
}
 \item{start.imp.sub}{
Starting value for the imputations of the outcome. When using binomial family, this is the value of the latent normal.
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

   \item{varY.start}{
    Starting value for varY, the residual variance of the substantive analysis model.  The default is the complete records estimate. When using family binomial, the variance is fixed to 1. 
  }
  
   \item{nburn}{
    Number of burn in iterations. Default is 1000.
  }
   \item{output}{
    When set to any value different from 1 (default), no output is shown on screen at the end of the process.
  }
  \item{out.iter}{
    When set to K, every K iterations a message "Iteration number N*K completed" is printed on screen. Default is 10. 
  }
   \item{family}{
    One of either "gaussian"" or "binomial". For binomial family, a probit link is assumed. 
  }
    
}

\value{
  A list is returned; this contains the final imputed dataset (finimp) and several 3-dimensional matrices, containing all the values drawn for each parameter at each iteration: these are fixed effect parameters of the covariates beta (collectbeta), level 1 covariance matrices (collectomega), fixed effect estimates of the substantive model and associated residual variances. If there are some categorical outcomes, a further output is included in the list, finimp.latnorm, containing the final state of the imputed dataset with the latent normal variables.
}

\examples{
   #First of all we load and attach the data:
  
  data(sldata)
  attach(sldata)
  
  #Then we make sure sex is a factor:
  
  sex<-factor(sex)
  
  # we define the data frame with all the variables 
  
  data<-data.frame(measure,age, sex)
  
  # And the formula of the substantive lm model 
  # sex as an outcome only because it is the only binary variable in the dataset...
  
  formula<-as.formula(sex~age+measure)
  
  #And finally we run the imputation function:
  
  imp<-jomo.glm.MCMCchain(formula,data, nburn=10)
  
  # Note we are using only 10 iterations to avoid time consuming examples,
  # which go against CRAN policies. In real applications we would use
  # much larger burn-ins (around 1000).
  
  # We can check, for example, the convergence of the first element of beta:
  
  plot(c(1:10),imp$collectbeta[1,1,1:10],type="l")
  
  
}
