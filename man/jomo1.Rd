\name{jomo1}
\alias{jomo1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JM Imputation of single level data
}
\description{
A wrapper function linking the 3 single level JM Imputation functions. The matrix of responses Y, must be a data.frame where continuous variables are numeric and binary/categorical variables are factors. 
}
\usage{
jomo1 (Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, 
      nburn=100, nbetween=100, nimp=5, output=1, out.iter=10)
}

\arguments{
  \item{Y}{
    A data.frame containing the outcomes of the imputation model. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
  }
 \item{X}{
    A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
  }
   \item{beta.start}{
    Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
 \item{l1cov.start}{
    Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
  }
   \item{l1cov.prior}{
    Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
  }
  \item{nburn}{
Number of burn in iterations. Default is 100.
}
  \item{nbetween}{
Number of iterations between two successive imputations. Default is 100.
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
}
\details{
This is just a wrapper function to link jomo1con, jomo1cat and jomo1mix. Format of the columns of Y is crucial in order for the function to be using the right sub-function. 
}
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 3-5, Wiley, ISBN: 978-0-470-74052-1.
}

\examples{

attach(sldata)

#We define all the inputs:

Y=data.frame(measure,age)
nburn=as.integer(200);
nbetween=as.integer(200);
nimp=as.integer(5);

# Then we run the function:

imp<-jomo1(Y,nburn=nburn,nbetween=nbetween,nimp=nimp)
}