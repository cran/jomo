\name{jomo1rancon.MCMCchain}
\alias{jomo1rancon.MCMCchain}

\title{
JM Imputation of clustered data with continuous variables only -
A tool to check convergence of the MCMC
}
\description{
This function is similar to jomo1rancon, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
jomo1rancon.MCMCchain(Y, X=NULL, Z=NULL, clus, beta.start=NULL, 
u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, 
l2cov.prior=NULL, start.imp=NULL, nburn=100, output=1, out.iter=10) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A data frame, or matrix, with responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
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
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. The default is a matrix of zeros.
}
  \item{u.start}{
A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
}
  \item{l1cov.start}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model. The default is the identity matrix.
}
  \item{l2cov.start}{
Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model times the number of random effects. The default is an identity matrix.
}
  \item{l1cov.prior}{
Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
}
  \item{l2cov.prior}{
Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
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
When set to K, every K iterations a message "Iteration number N*K completed" is printed on screen. Default is 10. 
}
}

\value{
A list with five elements is returned: the final imputed dataset (finimp) and four 3-dimensional matrices, containing all the values for beta (collectbeta), the random effects (collectu) and the level 1 (collectomega) and level 2 covariance matrices (collectcovu).
}

\examples{
#First of all we load and attach the data:
data(cldata)
attach(cldata)

#Then we define all the inputs:
Y<-data.frame(measure,age)
X<-data.frame(rep(1,1000),sex)
Z<-data.frame(rep(1,1000))
clus<-data.frame(city)
beta.start<-matrix(0,2,2)
u.start<-matrix(0,10,2)
l1cov.start<-diag(1,2)
l2cov.start<-diag(1,2)
l1cov.prior=diag(1,2);
nburn=as.integer(200);

l2cov.prior=diag(1,5);

#And finally we run the imputation function:
imp<-jomo1rancon.MCMCchain(Y,X,Z,clus,beta.start,u.start,l1cov.start, 
          l2cov.start,l1cov.prior,l2cov.prior,nburn=nburn)

#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")

#Or similarly we can check the convergence of any element of the level 2 covariance matrix:

plot(c(1:nburn),imp$collectcovu[1,2,1:nburn],type="l")

}
