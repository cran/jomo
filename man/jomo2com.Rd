\name{jomo2com}
\alias{jomo2com}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JM Imputation of 2-level data assuming a common level-1 covariance matrix across level-2 units.
}
\description{
Impute a 2-level dataset with mixed data types as outcome. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler where the covariance matrix is updated with a Metropolis-Hastings step. Fully observed categorical covariates may be considered as covariates as well, but they have to be included as dummy variables.
}
\usage{
jomo2com(Y.con=NULL, Y.cat=NULL, Y.numcat=NULL, Y2.con=NULL, Y2.cat=NULL,
Y2.numcat=NULL,X=NULL, X2=NULL, Z=NULL, clus, beta.start=NULL, l2.beta.start=NULL, 
u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, 
l2cov.prior=NULL, nburn=1000, nbetween=1000, nimp=5, output=1, out.iter=10) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.con}{
A data frame, or matrix, with level-1 continuous responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. 
}
  \item{Y.cat}{
A data frame, or matrix, with categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{Y.numcat}{
A vector with the number of categories in each categorical (or binary) variable.
}
  \item{Y2.con}{
A data frame, or matrix, with level-2 continuous responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. 
}
  \item{Y2.cat}{
A data frame, or matrix, with level-2 categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{Y2.numcat}{
A vector with the number of categories in each level-2 categorical (or binary) variable.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{X2}{
A data frame, or matrix, with level-2 covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{Z}{
A data frame, or matrix, for covariates associated to random effects in the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{clus}{
A data frame, or matrix, containing the cluster indicator for each observation. 
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
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
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
}
\details{
TThe Gibbs sampler algorithm used is described in detail in Chapter 9 of Carpenter and Kenward (2013). Regarding the choice of the priors, a flat prior is considered for beta and for the covariance matrix. A Metropolis Hastings step is implemented to update the covariance matrix, as described in the book. Binary or continuous covariates in the imputation model may be considered without any problem, but when considering a categorical covariate it has to be included with dummy variables (binary indicators) only.
}
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 9, Wiley, ISBN: 978-0-470-74052-1.
}

\examples{

Y<-tldata[,c("measure.a"), drop=FALSE]
Y2<-tldata[,c("big.city"), drop=FALSE]
clus<-tldata[,c("city")]

#now we run the imputation function. Note that we would typically use an higher 
#number of nburn iterations in real applications (at least 1000)

imp<-jomo2com(Y.con=Y, Y2.cat=Y2, Y2.numcat=2, clus=clus,nburn=10, nbetween=10, nimp=2)
  # Check help page for function jomo to see how to fit the model and 
  # combine estimates with Rubin's rules

}
