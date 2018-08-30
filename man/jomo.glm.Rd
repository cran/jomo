\name{jomo.glm}
\alias{jomo.glm}

\title{
  Joint Modelling Imputation Compatible with glm Model
}
\description{
  A function for substantive model compatible JM imputation, when the substantive model of interest is a simple generalized linear regression model. Interactions and polynomial functions of the covariates are allowed. Data must be passed as a data.frame where continuous variables are numeric and binary/categorical variables are factors.}
  
\usage{
 jomo.glm(formula, data, beta.start=NULL, l1cov.start=NULL, 
        l1cov.prior=NULL,nburn=1000, nbetween=1000, nimp=5, 
        output=1, out.iter=10, family="binomial") 
  
}

\arguments{

\item{formula}{
    an object of class formula: a symbolic description of the model to be fitted. It is possible to include in this formula interactions (through symbols '*' and '%') and polynomial terms (with the usual lm syntax, e.g. for a quadratic effect for variable x, 'I(x^2)'). 
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
    When set to K, every K iterations a dot is printed on screen. Default is 10. 
  }
   \item{family}{
    One of either "gaussian"" or "binomial". For binomial family, a probit link is assumed. 
  }
}
\details{
  This function allows for substantive model compatible imputation when the substantive model is a simple linear regression model. It can deal with interactions and polynomial terms through the usual lm syntax in the formula argument. 
  Format of the columns of data is crucial in order for the function to deal with binary/categorical covariates appropriately in the imputation algorithm.  }
\value{
  On screen, the posterior mean of the fixed effect estimates and of the residual variance are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
  Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Wiley, ISBN: 978-0-470-74052-1.
}

\examples{
   #First of all we attach the data (only for simplicity in these examples):
  
  attach(sldata)
  
  #Then we make sure sex is a factor:
  
  sex<-factor(sex)
  
  # we define the data frame with all the variables 
  
  data<-data.frame(measure,age, sex)
  
  # And the formula of the substantive lm model 
  # sex as an outcome only because it is the only binary variable in the dataset...
  
  formula<-as.formula(sex~age+measure)
  
  #And finally we run the imputation function:
  
  imp<-jomo.glm(formula,data, nburn=10, nbetween=10, nimp=2)
  
  # Note we are using only 10 iterations to avoid time consuming examples, 
  # which go against CRAN policies. In real applications we would use
  # much larger burn-ins (around 1000) and at least 5 imputations.
  
  
}