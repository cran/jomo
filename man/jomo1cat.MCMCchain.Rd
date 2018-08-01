\name{jomo1cat.MCMCchain}
\alias{jomo1cat.MCMCchain}

\title{
JM Imputation of single level data with categorical variables -
 A tool to check convergence of the MCMC
}
\description{
This function is similar to jomo1cat, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
jomo1cat.MCMCchain(Y.cat, Y.numcat, X=NULL, beta.start=NULL, 
l1cov.start=NULL, l1cov.prior=NULL, start.imp=NULL,
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
   \item{beta.start}{
    Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
  \item{l1cov.start}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
}
  \item{l1cov.prior}{
Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
}
 \item{start.imp}{
Starting value for the imputed dataset. n-level categorical variables are substituted by n-1 latent normals.
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
A list with four elements is returned: the final imputed dataset (finimp) and three 3-dimensional matrices, containing all the values drawn at each iteration for fixed effect parameters beta (collectbeta) and covariance matrix omega (collectomega). Finally, in finimp.latnorm, it is stored the final state of the imputed dataset with the latent normals in place of the categorical variables.

}

\examples{
attach(sldata)

# we define all the inputs:
# nimp, nburn and nbetween are smaller than they should. This is
#just because of CRAN policies on the examples.

Y.cat=data.frame(social)
Y.numcat=matrix(4,1,1)
X=data.frame(rep(1,300),sex)
beta.start<-matrix(0,2,3)
l1cov.start<-diag(1,3)
l1cov.prior=diag(1,3);
nburn=as.integer(100);

# Finally we run the sampler:

imp<-jomo1cat.MCMCchain(Y.cat,Y.numcat,X,beta.start,l1cov.start,l1cov.prior,nburn=nburn)

#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")


}
