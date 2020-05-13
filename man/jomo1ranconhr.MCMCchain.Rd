\name{jomo1ranconhr.MCMCchain}
\alias{jomo1ranconhr.MCMCchain}

\title{
JM Imputation of clustered data with continuous variables only with cluster-specific covariance matrices -
A tool to check convergence of the MCMC
}
\description{
This function is similar to jomo1ranconhr, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
jomo1ranconhr.MCMCchain(Y, X=NULL, Z=NULL, clus, 
beta.start=NULL, u.start=NULL, l1cov.start=NULL, 
l2cov.start=NULL, l1cov.prior=NULL, l2cov.prior=NULL,start.imp=NULL,  
nburn=1000, a=(ncol(Y)+50),a.prior=NULL, meth="random", output=1, out.iter=10) 

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
  \item{l1cov.start}{ in column. Dimension of each square matrix is equal to the number of outcomes in the imputation model. The default is the identity matrix for each cluster.
}
  \item{l2cov.start}{
Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model times the number of random effects. The default is an identity matrix.
}
  \item{l1cov.prior}{
Scale matrix for the inverse-Wishart prior for the covariance matrices. The default is the identity matrix.
}
  \item{l2cov.prior}{
Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
}
\item{start.imp}{
Starting value for the imputed dataset. 
}
   \item{nburn}{
Number of iterations. Default is 1000.
}

  \item{a}{
Starting value for the degrees of freedom of the inverse Wishart distribution of the cluster-specific covariance matrices. Default is 50+D, with D being the dimension of the covariance matrices.
}
  \item{a.prior}{
Hyperparameter (Degrees of freedom) of the chi square prior distribution for the degrees of freedom of the inverse Wishart distribution for the cluster-specific covariance matrices. Default is D, with D being the dimension of the covariance matrices.
}

\item{meth}{
This can be set to "Fixed" or "Random". In the first case the function will consider fixed study-specific covariance matrices, in the second, random study-specific distributed according to an inverse-Wishart distribution.
}
 \item{output}{
When set to any value different from 1 (default), no output is shown on screen at the end of the process.
}
\item{out.iter}{
When set to K, every K iterations a dot is printed on screen. Default is 10. 
}
}

\value{
A list with five elements is returned: the final imputed dataset (finimp) and four 3-dimensional matrices, containing all the values for beta (collectbeta), the random effects (collectu) and the level 1 (collectomega) and level 2 covariance matrices (collectcovu).
}

\examples{



# we define the inputs
# nburn is smaller than needed. This is
#just because of CRAN policies on the examples.

Y<-cldata[,c("measure","age")]
clus<-cldata[,c("city")]
X=data.frame(rep(1,1000),cldata[,c("sex")])
colnames(X)<-c("const", "sex")
Z<-data.frame(rep(1,1000))
nburn=as.integer(200);
a=3

# Finally we run either the model with fixed or random cluster-specific cov. matrices:

imp<-jomo1ranconhr.MCMCchain(Y,X,Z,clus,nburn=nburn,meth="random")
          
#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")

#Or similarly we can check the convergence of any element of the level 2 cov. matrix:

plot(c(1:nburn),imp$collectcovu[1,2,1:nburn],type="l")
}
