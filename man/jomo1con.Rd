\name{jomo1con}
\alias{jomo1con}
 
\title{
JM Imputation of single level data with continuous variables only
}
\description{
Impute a single level dataset with continuous outcomes only. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler. Categorical covariates may be considered, but they have to be included with dummy variables.
}
\usage{
jomo1con(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL, 
nburn=100, nbetween=100, nimp=5, output=1, out.iter=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A data frame, or matrix, with responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{beta.start}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. The default is a matrix of zeros.
}
  \item{l1cov.start}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model. The default is the identity matrix.
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
The Gibbs sampler algorithm used is described in detail in Chapter 3 of Carpenter and Kenward (2013). Regarding the choice of the priors, a flat prior is considered for beta, while an inverse-Wishart prior is given to the covariance matrix, with p-1 degrees of freedom, aka the minimum possible, to guarantee the greatest uncertainty. Binary or continuous covariates in the imputation model may be considered without any problem, but when considering a categorical covariate it has to be included through dummy variables (binary indicators) only. 
}
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 3, Wiley, ISBN: 978-0-470-74052-1.
}

\examples{

#We define all the inputs:

Y=sldata[,c("measure", "age")]
X=data.frame(rep(1,300),sldata[,c("sex")])
colnames(X)<-c("const", "sex")
beta.start<-matrix(0,2,2)
l1cov.start<-diag(1,2)
l1cov.prior=diag(1,2);
nburn=as.integer(200);
nbetween=as.integer(200);
nimp=as.integer(5);

# Then we run he function:

imp<-jomo1con(Y,X,beta.start,l1cov.start,l1cov.prior,nburn,nbetween,nimp)

cat("Original value was missing(",imp[1,1],"), imputed value:", imp[301,1])

  # Check help page for function jomo to see how to fit the model and 
  # combine estimates with Rubin's rules


}
