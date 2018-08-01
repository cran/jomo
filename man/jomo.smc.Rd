\name{jomo.smc}
\alias{jomo.smc}

\title{
  Joint Modelling Substantive Model Compatible Imputation 
}
\description{
  A wrapper function for all the substantive model compatible JM imputation functions. The substantive model of interest is either lm, glm, lmer, glmer or coxph. Interactions and polynomial functions of the covariates are allowed. Data must be passed as a data.frame where continuous variables are numeric and binary/categorical variables are factors.}

\usage{
  jomo.smc(formula, data, level=rep(1,ncol(data)), beta.start=NULL,
  l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL,
  l1cov.prior=NULL, l2cov.prior=NULL, a.start=NULL, a.prior=NULL, 
  nburn=1000, nbetween=1000, nimp=5, meth="common", family="binomial",
  output=1, out.iter=10, model) 
    
}

\arguments{
  
  \item{formula}{
    an object of class formula: a symbolic description of the model to be fitted. It is possible to include in this formula interactions (through symbols '*' and '%') and polynomial terms (with the usual lm syntax, e.g. for a quadratic effect for variable x, 'I(x^2)'). Random effects follow the usual lmer syntax; however, it is possible only to allow for one grouping variable and all the variables with random effects must be included within the same brackets (es: (1+X|clus) is correct, while (1|clus)+(X|clus) is NOT).
  }
  
  \item{data}{
    A data.frame containing all the variables to include in the imputation model. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
  }
  \item{level}{
  If the dataset is multilevel, this must be a vector indicating whether each variable is either a level 1 or a level 2 variable. The value assigned to the cluster indicator is irrelevant.
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

  \item{nburn}{
    Number of burn in iterations. Default is 1000.
  }
  \item{nbetween}{
    Number of iterations between two successive imputations. Default is 1000.
  }
  \item{nimp}{
    Number of Imputations. Default is 5.
  }
  \item{output}{
    When set to any value different from 1 (default), no output is shown on screen at the end of the process.
  }
  \item{out.iter}{
    When set to K, every K iterations a message "Iteration number N*K completed" is printed on screen. Default is 10. 
  }
   \item{model}{
    The type of model we want to impute compatibly with. It can currently be one of lm, glm (binomial), coxph, lmer or glmer (binomial).
  }
     \item{family}{
    One of either "gaussian"" or "binomial". For binomial family, a probit link is assumed. 
  }
  
}
\details{
  This function allows for substantive model compatible imputation. It can deal with interactions and polynomial terms through the usual lmer syntax in the formula argument. 
  Format of the columns of data is crucial in order for the function to deal with binary/categorical covariates appropriately in the imputation algorithm.  }
\value{
  On screen, the posterior mean of the fixed effect estimates and of the residual variance are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
  Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Wiley, ISBN: 978-0-470-74052-1.
}

\examples{

  attach(cldata)
  
  # make sure sex is a factor:
  
  sex<-factor(sex)
  
  # we define the data frame with all the variables 
  
  data<-data.frame(measure,age, sex, city)
  mylevel<-c(1,1,1,1)
  
  # And the formula of the substantive lm model
  
  formula<-as.formula(measure~sex+age+I(age^2)+(1|city))
  
  #And finally we run the imputation function:
  
  imp<-jomo.smc(formula,data, level=mylevel, nburn=10, nbetween=10, model="lmer")
  
  # Note we are using only 10 iterations to avoid time consuming examples, 
  # which go against CRAN policies. 
  # If we were interested in a model with interactions:
  
  # formula2<-as.formula(measure~sex*age+(1|city))
  # imp2<-jomo.smc(formula2,data, level=mylevel, nburn=10, nbetween=10, model="lmer")
  
  # The analysis and combination steps are as for all the other functions
  # (see e.g. help file for jomo1con)
  
  
  
}