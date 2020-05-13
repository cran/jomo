\name{jomo.coxph}
\alias{jomo.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Joint Modelling Imputation Compatible with Cox Proportional Hazards Model
}
\description{
  A function for substantive model compatible JM imputation, when the substantive model of interest is a Cox Proportional Hazards Model. Interactions and polynomial functions of the covariates are allowed. Data must be passed as a data.frame where continuous variables are numeric and binary/categorical variables are factors.}

\usage{
  jomo.coxph(formula, data,  beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, 
          nburn=1000, nbetween=1000, nimp=5, output=1, out.iter=10) 
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
    When set to 0, no output is shown on screen at the end of the process. When set to 1, only the parameter estimates related to the substantive model are shown (default). When set to 2, all parameter estimates (posterior means) are displayed.
  }
  \item{out.iter}{
    When set to K, every K iterations a dot is printed on screen. Default is 10. 
  }

  
}
\details{ This function allows for substantive model compatible imputation when the substantive model is a Cox PH model. It can deal with interactions and polynomial terms through the usual lm syntax in the formula argument. 
  Format of the columns of data is crucial in order for the function to deal with binary/categorical covariates appropriately in the imputation algorithm.  }

\value{
  On screen, the posterior mean of the fixed effect estimates and of the residual variance are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}

\examples{
#define substantive model
formula<-as.formula(Surv(time, status) ~ measure + sex + I(measure^2))

#Run imputation
if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  imp<-jomo.coxph(formula,surdata, nburn = 100, nbetween = 100, nimp=5)
}
  # Check help page for function jomo to see how to fit the model and 
  # combine estimates with Rubin's rules

}

