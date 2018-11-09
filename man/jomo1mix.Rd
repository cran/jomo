\name{jomo1mix}
\alias{jomo1mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JM Imputation of single level data with mixed variable types 
}
\description{
Impute a single level dataset with mixed data types as outcome. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler where the covariance matrix is updated  with a Metropolis-Hastings step. Fully observed categorical variables may be considered as covariates as well, but they have to be included as dummy variables.
}
\usage{
jomo1mix(Y.con, Y.cat, Y.numcat, X=NULL, beta.start=NULL, l1cov.start=NULL, 
l1cov.prior=NULL, nburn=100, nbetween=100, nimp=5, output=1,out.iter=10)
 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.con}{
A data frame, or matrix, with continuous responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA. If no continuous outcomes are present in the model, jomo1cat should be used instead.
}
  \item{Y.cat}{
A data frame, or matrix, with categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{Y.numcat}{
A vector with the number of categories in each categorical (or binary) variable.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{beta.start}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we define n-1 latent normals. The default is a matrix of zeros.
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
Regarding the choice of the priors, a flat prior is considered for beta and for the covariance matrix. A Metropolis Hastings step is implemented to update the covariance matrix, as described in the book. Binary or continuous covariates in the imputation model may be considered without any problem, but when considering a categorical covariate it has to be included with dummy variables (binary indicators) only. 
}
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 5, Wiley, ISBN: 978-0-470-74052-1.
}

\examples{

#Then, we define all the inputs:
# nburn is smaller than needed. This is
#just because of CRAN policies on the examples.

Y.con=sldata[,c("measure","age")]
Y.cat=sldata[,c("social"), drop=FALSE]
Y.numcat=matrix(4,1,1)
X=data.frame(rep(1,300),sldata[,c("sex")])
colnames(X)<-c("const", "sex")
beta.start<-matrix(0,2,5)
l1cov.start<-diag(1,5)
l1cov.prior=diag(1,5);
nburn=as.integer(100);
nbetween=as.integer(100);
nimp=as.integer(5);

#Then we run the sampler:

imp<-jomo1mix(Y.con,Y.cat,Y.numcat,X,beta.start,l1cov.start,
      l1cov.prior,nburn,nbetween,nimp)

cat("Original value was missing(",imp[1,1],"), imputed value:", imp[301,1])

  # Check help page for function jomo to see how to fit the model and 
  # combine estimates with Rubin's rules


}
