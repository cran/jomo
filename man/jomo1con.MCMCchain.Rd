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
jomo1con.MCMCchain(Y, X=matrix(1,nrow(Y),1), 
betap=matrix(0,ncol(X),ncol(Y)), covp=diag(1,ncol(Y)), 
Sp=diag(1,ncol(Y)), nburn=100, output=1, out.iter=10) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A data frame, or matrix, with responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{betap}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. The default is a matrix of zeros.
}
  \item{covp}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model. The default is the identity matrix.
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
# Using sldata:

data(sldata)
attach(sldata)

#We define all the inputs:

Y=data.frame(measure,age)
X=data.frame(rep(1,300),sex)
betap<-matrix(0,2,2)
covp<-diag(1,2)
Sp=diag(1,2);
nburn=as.integer(200);

# Then we run he function:

imp<-jomo1con.MCMCchain(Y,X,betap,covp,Sp,nburn)

#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")

#Or similarly we can check the convergence of any element of omega:

plot(c(1:nburn),imp$collectomega[1,2,1:nburn],type="l")
}
