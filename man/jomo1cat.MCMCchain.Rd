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
jomo1cat.MCMCchain(Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), 
betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat)))), covp=diag(1,ncol(betap)), 
Sp=diag(1,ncol(betap)), nburn=100, output=1, out.iter=10) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{Y_cat}{
A data frame, or matrix, with categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{Y_numcat}{
A vector with the number of categories in each categorical (or binary) variable.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{betap}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we define n-1 latent normals. The default is a matrix of zeros.
}
  \item{covp}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
}
  \item{Sp}{
Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
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
A list with three elements is returned: the final imputed dataset (finimp) and three 3-dimensional matrices, containing all the values for beta (collectbeta) and omega (collectomega).

}

\examples{
#First of all we load and attach sldata
library(jomo)
data(sldata)
attach(sldata)

#Then, we define all the inputs:
# nimp, nburn and nbetween are smaller than they should. This is
#just because of CRAN policies on the examples.

Y_cat=data.frame(social)
Y_numcat=matrix(4,1,1)
X=data.frame(rep(1,300),sex)
betap<-matrix(0,2,3)
covp<-diag(1,3)
Sp=diag(1,3);
nburn=as.integer(100);

# Finally we run the sampler:

imp<-jomo1cat.MCMCchain(Y_cat,Y_numcat,X,betap,covp,Sp,nburn)

#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")


}
