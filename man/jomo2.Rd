\name{jomo2}
\alias{jomo2}

\title{
JM Imputation of 2-level data
}
\description{
A wrapper function linking the 2-level JM Imputation functions. The matrices of responses Y and Y2, must be data.frames where continuous variables are numeric and binary/categorical variables are factors.}
\usage{
jomo2(Y, Y2, X=NULL, X2=NULL, Z=NULL,clus, beta.start=NULL, l2.beta.start=NULL,
u.start=NULL, l1cov.start=NULL, l2cov.start=NULL, l1cov.prior=NULL, 
l2cov.prior=NULL, nburn=1000, nbetween=1000, nimp=5, a=NULL, a.prior=NULL,
meth="common", output=1, out.iter=10) 

}

\arguments{
  \item{Y}{
A data.frame with the level-1 outcomes of the imputation model, where columns related to continuous variables are numeric and columns related to binary/categorical variables are factors. 
}
   \item{Y2}{
    A data.frame containing the level-2 outcomes of the imputation model, i.e. the partially observed level-2 variables. Columns related to continuous variables have to be numeric and columns related to binary/categorical variables have to be factors. 
  }
  \item{X}{
    A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
  }
  \item{X2}{
    A data frame, or matrix, with level-2 covariates of the joint imputation model. Rows correspond to different level-1 observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
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
  \item{nburn}{
Number of burn in iterations. Default is 1000.
}
  \item{nbetween}{
Number of iterations between two successive imputations. Default is 1000.
}
  \item{nimp}{
Number of Imputations. Default is 5.
}
  \item{a}{
Starting value for the degrees of freedom of the inverse Wishart distribution of the cluster-specific covariance matrices. Default is 50+D, with D being the dimension of the covariance matrices. This is used only when option meth is set to "random".
}
  \item{a.prior}{
Hyperparameter (Degrees of freedom) of the chi square prior distribution for the degrees of freedom of the inverse Wishart distribution for the cluster-specific covariance matrices. Default is D, with D being the dimension of the covariance matrices.
}

\item{meth}{
Method used to deal with level 1 covariance matrix. When set to "common", a common matrix across clusters is used (function jomo2com). When set to "fixed", fixed study-specific matrices are considered (jomo2hr with option meth="fixed"). Finally, when set to "random", random study-specific matrices are considered (jomo2hr with option meth="random")
}

 \item{output}{
When set to any value different from 1 (default), no output is shown on screen at the end of the process.
}
\item{out.iter}{
When set to K, every K iterations a dot is printed on screen. Default is 10. 
}

}
\details{
This is just a wrapper function to link jomo1rancon, jomo1rancat and jomo1ranmix and the respective "hr" (heterogeneity in covariance matrices) versions. Format of the columns of Y is crucial in order for the function to be using the right sub-function. }
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 9, Wiley, ISBN: 978-0-470-74052-1.
}

\examples{

Y<-tldata[,c("measure.a"), drop=FALSE]
Y2<-tldata[,c("big.city"), drop=FALSE]
clus<-tldata[,c("city")]
nburn=3
nbetween=3
nimp=2

#now we run the imputation function. Note that we would typically use an higher 
#number of nburn iterations in real applications (at least 1000)

imp<-jomo2(Y=Y, Y2=Y2, clus=clus,nburn=nburn, nbetween=nbetween, nimp=nimp)

  # Check help page for function jomo to see how to fit the model and 
  # combine estimates with Rubin's rules


}