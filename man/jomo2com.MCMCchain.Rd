\name{jomo2com.MCMCchain}
\alias{jomo2com.MCMCchain}

\title{
JM Imputation of 2-level data assuming a common level-1 covariance matrix across level-2 units -
A tool to check convergence of the MCMC
}
\description{
This function is similar to jomo2com, but it returns the values of all the parameters in the model at each step of the MCMC instead of the imputations. It is useful to check the convergence of the MCMC sampler.
}
\usage{
jomo2com.MCMCchain(Y.con=NULL, Y.cat=NULL, Y.numcat=NULL, Y2.con=NULL, 
Y2.cat=NULL, Y2.numcat=NULL, X=NULL, X2=NULL, Z=NULL, clus, beta.start=NULL,
l2.beta.start=NULL, u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, 
l1cov.prior=NULL, l2cov.prior=NULL, start.imp=NULL, l2.start.imp=NULL, nburn=100, 
output=1, out.iter=10) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.con}{
A data frame, or matrix, with level-1 continuous responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. 
}
  \item{Y.cat}{
A data frame, or matrix, with categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{Y.numcat}{
A vector with the number of categories in each categorical (or binary) variable.
}
  \item{Y2.con}{
A data frame, or matrix, with level-2 continuous responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. 
}
  \item{Y2.cat}{
A data frame, or matrix, with level-2 categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA.
}
  \item{Y2.numcat}{
A vector with the number of categories in each level-2 categorical (or binary) variable.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{X2}{
A data frame, or matrix, with level-2 covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{Z}{
A data frame, or matrix, for covariates associated to random effects in the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{clus}{
A data frame, or matrix, containing the cluster indicator for each observation. 
}
   \item{beta.start}{
    Starting value for beta, the vector(s) of level-1 fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
  \item{l2.beta.start}{
    Starting value for beta2, the vector(s) of level-2 fixed effects. Rows index different covariates and columns index different level-2 outcomes. For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals. The default is a matrix of zeros.
  }
  \item{u.start}{
A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
}
  \item{l1cov.start}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
}
  \item{l2cov.start}{
    Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model times the number of random effects plus the number of level-2 outcomes. The default is an identity matrix.
  }
  \item{l1cov.prior}{
Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
}
  \item{l2cov.prior}{
Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
}
 \item{start.imp}{
Starting value for the imputed dataset. n-level categorical variables are substituted by n-1 latent normals.
}
  \item{l2.start.imp}{
Starting value for the level-2 imputed variables. n-level categorical variables are substituted by n-1 latent normals.
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
 A list is returned; this contains the final imputed dataset (finimp) and several 3-dimensional matrices, containing all the values drawn for each parameter at each iteration: these are, potentially, fixed effect parameters beta (collectbeta), random effects (collectu), level 1 (collectomega) and level 2 covariance matrices (collectcovu) and level-2 fixed effect parameters. If there are some categorical outcomes, a further output is included in the list, finimp_latnorm, containing the final state of the imputed dataset with the latent normal variables.
}

\examples{
 #First of all we load and attach the data:
  
data(tldata)
attach(tldata)
nburn=20

#And then we run the imputation function. Note that we would typically use an higher 
#number of nburn iterations in real applications (at least 100)

imp<-jomo2com.MCMCchain(Y.con=data.frame(measure.a), 
          Y2.cat=data.frame(big.city), Y2.numcat=c(2), clus=city,nburn=nburn)
          
#We can check the convergence of the first element of beta:

plot(c(1:nburn),imp$collectbeta[1,1,1:nburn],type="l")
}
