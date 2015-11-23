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
jomo1ranconhr.MCMCchain(Y, X=matrix(1,nrow(Y),1), Z=matrix(1,nrow(Y),1), clus, 
betap=matrix(0,ncol(X),ncol(Y)), up=matrix(0,nrow(unique(clus)),ncol(Z)*ncol(Y)), 
covp=matrix(diag(1,ncol(Y)),nrow(unique(clus))*ncol(Y),ncol(Y),2), 
covu=diag(1,ncol(Y)*ncol(Z)), Sp=diag(1,ncol(Y)), Sup=diag(1,ncol(Y)*ncol(Z)), 
nburn=100, a=ncol(Y),meth="random", output=1, out.iter=10) 

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
  \item{betap}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. The default is a matrix of zeros.
}
  \item{up}{
A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
}
  \item{covp}{
Starting value for the covariance matrices, pulled one above the other in column. Dimension of each square matrix is equal to the number of outcomes in the imputation model. The default is the identity matrix for each cluster.
}
  \item{covu}{
Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model times the number of random effects. The default is an identity matrix.
}
  \item{Sp}{
Scale matrix for the inverse-Wishart prior for the covariance matrices. The default is the identity matrix.
}
  \item{Sup}{
Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
}
   \item{nburn}{
Number of iterations. Default is 100.
}

  \item{a}{
Starting value for the degrees of freedom of the wishart distribution from which all of the covariance matrices are drawn in case method="random" option in specified. Default is the minimum possible, i.e. the dimension of the covariance matrices.
}
\item{meth}{
This can be set to "Fixed" or "Random". In the first case the function will consider fixed study-specific covariance matrices, in the second, random study-specific distributed according to an inverse-Wishart distribution.
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

#Then we define the inputs
# nimp, nburn and nbetween are smaller than they should. This is
#just because of CRAN policies on the examples.

Y<-data.frame(measure,age)
X<-data.frame(rep(1,1000),sex)
Z<-data.frame(rep(1,1000))
clus<-data.frame(city)
betap<-matrix(0,2,2)
up<-matrix(0,10,2)
covp<-matrix(diag(1,2),20,2,2)
covu<-diag(1,2)
Sp=diag(1,2);
nburn=as.integer(50);

Sup=diag(1,5);
a=3

# And finally we can run either the model with fixed or random cluster-specific covariance matrices:

imp<-jomo1ranconhr.MCMCchain(Y,X,Z,clus,betap,up,covp, covu,Sp,Sup,nburn,meth="random")
#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")

#Or similarly we can check the convergence of any element of th elevel 2 covariance matrix:

plot(c(1:nburn),imp$collectcovu[1,2,1:nburn],type="l")
}
