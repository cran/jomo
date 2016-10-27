\name{jomo1.MCMCchain}
\alias{jomo1.MCMCchain}

\title{
JM Imputation of single level data - A tool to check convergence of the MCMC
}
\description{
This function is similar to jomo1, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
jomo1.MCMCchain(Y, X=NULL, beta.start=NULL, l1cov.start=NULL, l1cov.prior=NULL,
start.imp=NULL, nburn=100, output=1, out.iter=10)
}

\arguments{
  \item{Y}{
    A data.frame containing the outcomes of the imputation model. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
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
A list with three elements is returned: the final imputed dataset (finimp) and three 3-dimensional matrices, containing all the values for beta (collectbeta) and omega (collectomega). If there are some categorical outcomes, a further output is included in the list, finimp.latnorm, containing the final state of the imputed dataset with the latent normal variables.
}

\examples{
# Using sldata:

data(sldata)
attach(sldata)

#We define all the inputs:

Y=data.frame(measure,age)
nburn=as.integer(200);

# Then we run the function:

imp<-jomo1.MCMCchain(Y,nburn=nburn)

#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")

#Or similarly we can check the convergence of any element of omega:

plot(c(1:nburn),imp$collectomega[1,2,1:nburn],type="l")


}
