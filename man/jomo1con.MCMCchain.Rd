\name{jomo1con.MCMCchain}
\alias{jomo1con.MCMCchain}

\title{
JM Imputation of single level data with continuous variables only -
A tool to check convergence of the MCMC
}
\description{
This function is similar to jomo1con, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
jomo1con.MCMCchain(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, 
l1cov.prior=NULL, start.imp=NULL, nburn=100, output=1, out.iter=10) 

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
 \item{start.imp}{
Starting value for the imputed dataset. 
}
  \item{nburn}{
Number of iterations. Default is 100.
}

 \item{output}{
When set to any value different from 1 (default), no output is shown on screen at the end of the process.
}
\item{out.iter}{
When set to K, every K iterations a dot is printed on screen. Default is 10. 
}
}

\value{
A list with three elements is returned: the final imputed dataset (finimp) and three 3-dimensional matrices, containing all the values for the fixed effect parameters beta (collectbeta) and the covariance matrix omega (collectomega).
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

# Then we run he function:

imp<-jomo1con.MCMCchain(Y,X,beta.start,l1cov.start,l1cov.prior,nburn=nburn)

#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")

#Or similarly we can check the convergence of any element of omega:

plot(c(1:nburn),imp$collectomega[1,2,1:nburn],type="l")
}
