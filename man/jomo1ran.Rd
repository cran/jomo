\name{jomo1ran}
\alias{jomo1ran}

\title{
JM Imputation of 2-level data
}
\description{
A wrapper function linking the six 2-level JM Imputation functions. The matrix of responses Y, must be a data.frame where continuous variables are numeric and binary/categorical variables are factors.}
\usage{
jomo1ran(Y, X=NULL, Z=NULL,clus, 
      betap=NULL, up=NULL, covp=NULL, covu=NULL, Sp=NULL, Sup=NULL, 
      nburn=500, nbetween=100, nimp=5, a=NULL, meth="common") 

}

\arguments{
   \item{Y}{
a data.frame where columns related to continuous variables are numeric and columns related to binary/categorical variables are factors. 
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{Z}{
A data frame, or matrix, for covariates associated to random effects in the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{clus}{
A data frame, or matrix, containing the cluster indicator for each observation. Cluster needs to be labeled with an integer number ranging from 0 to nclus-1.
}
  \item{betap}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we define n-1 latent normals. The default is a matrix of zeros.
}
  \item{up}{
A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
}
  \item{covp}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix.
}
  \item{covu}{
Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model times the number of random effects. The default is an identity matrix.
}
  \item{Sp}{
Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
}
  \item{Sup}{
Scale matrix for the inverse-Wishart prior for the level 2 covariance matrix. The default is the identity matrix.
}
  \item{nburn}{
Number of burn in iterations. Default is 100.
}
  \item{nbetween}{
Number of iterations between two successive imputations. Default is 100.
}
  \item{nimp}{
Number of Imputations. Default is 5.
}
\item{a}{
Starting value for the degrees of freedom of the inverse Wishart distribution from which all of the covariance matrices are drawn in case method="random" option in specified. Default is the minimum possible, i.e. the dimension of the covariance matrices.
}

\item{meth}{
Method used to deal with level 1 covariance matrix. When set to "common", a common matrix across clusters is used (functions jomo1rancon, jomo1rancat and jomo1ranmix). When set to "fixed", fixed study-specific matrices are considered (jomo1ranconhr, jomo1rancathr and jomo1ranmixhr with coption meth="fixed"). Finally, when set to "random", random study-specific matrices are considered (jomo1ranconhr, jomo1rancathr and jomo1ranmixhr with coption meth="random")
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
  #First of all we load and atach the data:
  
data(mldata)
attach(mldata)

#Then we define all the inputs:

Y<-data.frame(measure,age)
clus<-data.frame(city)
nburn=as.integer(200);
nbetween=as.integer(200);
nimp=as.integer(5);


#And finally we run the imputation function:
imp<-jomo1ran(Y,clus=clus,nburn=nburn,nbetween=nbetween,nimp=nimp)

#we could even run it with fixed or random cluster-specific covariance matrices:

#imp<-jomo1ran(Y,clus=clus,nburn=nburn,nbetween=nbetween,nimp=nimp, meth="fixed")
#or:
#imp<-jomo1ran(Y,clus=clus,nburn=nburn,nbetween=nbetween,nimp=nimp, meth="random")

}