\name{jomo1rancat.MCMCchain}
\alias{jomo1rancat.MCMCchain}

\title{
JM Imputation of clustered data with categorical variables - 
A tool to check convergence of the MCMC
}
\description{
This function is similar to jomo1rancat, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
jomo1rancat.MCMCchain(Y.cat, Y.numcat, X=matrix(1,nrow(Y.cat),1), 
Z=matrix(1,nrow(Y.cat),1),clus, 
beta.start=matrix(0,ncol(X),((sum(Y.numcat)-length(Y.numcat)))), 
u.start=matrix(0,nrow(unique(clus)),
    ncol(Z)*((sum(Y.numcat)-length(Y.numcat)))), 
l1cov.start=diag(1,ncol(beta.start)), l2cov.start=diag(1,ncol(u.start)), 
l1cov.prior=diag(1,ncol(l1cov.start)), l2cov.prior=diag(1,ncol(l2cov.start)), 
nburn=100, output=1, out.iter=10) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.cat}{
A data frame, or matrix, with categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{Y.numcat}{
A vector with the number of categories in each categorical (or binary) variable.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{Z}{
A data frame, or matrix, for covariates associated to random effects in the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{clus}{
A data frame, or matrix, containing the cluster indicator for each observation. 
}
  \item{beta.start}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we define n-1 latent normals. The default is a matrix of zeros.
}
  \item{u.start}{
A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
}
  \item{l1cov.start}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
}
  \item{l2cov.start}{
Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model times the number of random effects. The default is an identity matrix.
}
  \item{l1cov.prior}{
Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
}
  \item{l2cov.prior}{
Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
}
  \item{nburn}{
Number of burn in iterations. Default is 100.
}
 \item{output}{
When set to any value different from 1 (default), no output is shown on screen at the end of the process.
}
\item{out.iter}{
When set to K, every K iterations a message "Iteration number N*K completed" is printed on screen. Default is 10. 
}
}
\value{
A list with five elements is returned: the final imputed dataset (finimp) and four 3-dimensional matrices, containing all the values for beta (collectbeta), the random effects (collectu) and the level 1 (collectomega) and level 2 covariance matrices (collectcovu).
}

\examples{
#First of all we load and attach the data:
data(mldata)
attach(mldata)

#Then we define all the inputs:
# nimp, nburn and nbetween are smaller than they should. This is
#just because of CRAN policies on the examples.

Y.cat=data.frame(social)
Y.numcat=matrix(4,1,1)
X=data.frame(rep(1,1000),sex)
Z<-data.frame(rep(1,1000))
clus<-data.frame(city)
beta.start<-matrix(0,2,3)
u.start<-matrix(0,10,3)
l1cov.start<-diag(1,3)
l2cov.start<-diag(1,3)
l1cov.prior=diag(1,3);
l2cov.prior=diag(1,3);
nburn=as.integer(100);

#And finally we run the imputation function:

imp<-jomo1rancat.MCMCchain(Y.cat, Y.numcat, X,Z,clus,beta.start,u.start,l1cov.start, 
l2cov.start,l1cov.prior,l2cov.prior,nburn)
#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")

#Or similarly we can check the convergence of any element of th elevel 2 covariance matrix:

plot(c(1:nburn),imp$collectcovu[1,2,1:nburn],type="l")
}