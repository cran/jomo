\name{jomo1rancathr}
\alias{jomo1rancathr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JM Imputation of clustered data with categorical variables with cluster-specific covariance matrices
}
\description{
Impute a clustered dataset with categorical variables as outcome. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler where a different covariance matrix is sampled within each cluster from the same inverse Wishart distribution. Fully observed categorical covariates may be considered as covariates as well, but they have to be included as dummy variables.
}
\usage{
jomo1rancathr( Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), 
Z=matrix(1,nrow(Y_cat),1), clus,
betap=matrix(0,ncol(X),((sum(Y_numcat)-length(Y_numcat)))), 
up=matrix(0,length(unique(clus)),ncol(Z)*((sum(Y_numcat)-length(Y_numcat)))), 
covp=matrix(diag(1,ncol(betap)),ncol(betap)*length(unique(clus)),ncol(betap),2), 
covu=diag(1,ncol(up)), Sp=diag(1,ncol(betap)), Sup=diag(1,ncol(up)), 
nburn=100, nbetween=100, nimp=5,a=ncol(betap), rngflag=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Y_cat}{
A data frame, or matrix, with categorical (or binary) responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Categories must be integer numbers from 1 to N. Missing values are coded as NA.
}
  \item{Y_numcat}{
A vector with the number of categories in each categorical (or binary) variable.
}
  \item{X}{
A data frame, or matrix, with covariates of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{Z}{
A data frame, or matrix, for covariates associated to random effects in the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are not allowed in these variables. In case we want an intercept, a column of 1 is needed. The default is a column of 1.
}
  \item{clus}{
A vector containing the cluster indicator for each observation. Cluster needs to be labeled with an integer number ranging from 0 to nclus-1.
}
  \item{betap}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. For each n-category variable we define n-1 latent normals. The default is a matrix of zeros.
}
  \item{up}{
A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
}
  \item{covp}{
Starting value for the covariance matrices, pulled one above the other in column. Dimension of each square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model. The default is the identity matrix for each cluster.
}
  \item{covu}{
Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes (continuous plus latent normals) in the imputation model times the number of random effects. The default is an identity matrix.
}
   \item{Sp}{
Scale matrix for the inverse-Wishart prior for the covariance matrices. The default is the identity matrix.
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
Starting value for the degrees of freedom of the wishart distribution from which all of the covariance matrices are drawn. Default is the minimum possible, i.e. the dimension of the covariance matrices.
}
  \item{rngflag}{
A flag for the random number generator. If set to zero (default) rnglib is used to initialize, from library rnglib.
Otherwise R random number generator is used. Suggestion is to use rnglib becaus eof its speed, but R rng when looking for reproducible results.
}
}
\details{
The Gibbs sampler algorithm used is obtained is a mixture of the ones described in chapter 5 and 9 of Carpenter and Kenward (2013). We update the covariance matrices drawing from the inverse Wishart distribution as we usually do when only continuous data are present in the model, and afterwards we constrain variances of the latent normals to 1 and covariances among latent normals related to the same categorical variable to 0. Then positive definitiveness of the matrix is checked and, if not met, another matrix is drawn. We also update values of a and A, degrees of freedom and scale matrix of the inverse Wishart distribution from which all the covariance matrices are sampled. Regarding the choice of the priors, a flat prior is considered for beta, while an inverse-Wishart prior is given to the covariance matrix, with p-1 degrees of freedom, aka the minimum possible, to guarantee the greatest uncertainty.   Binary or continuous covariates in the imputation model may be considered without any problem, but when considering a categorical covariate it has to be included with dummy variables (binary indicators) only. 
}
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 9, Wiley, ISBN: 978-0-470-74052-1.

Yucel R.M., (2011), Random-covariances and mixed-effects models for imputing multivariate multilevel continuous data, Statistical Modelling, 11 (4), 351-370, DOI: 10.1177/1471082X100110040.

}

\examples{

#First of all we load and attach the data:

data(mldata)
attach(mldata)

#Then we define the inputs
# nimp, nburn and nbetween are smaller than they should. This is
#just because of CRAN policies on the examples.

Y_cat=data.frame(social)
Y_numcat=matrix(4,1,1)
X=data.frame(rep(1,1000),sex)
Z<-data.frame(rep(1,1000))
clus<-data.frame(city)
betap<-matrix(0,2,3)
up<-matrix(0,10,3)
covp<-matrix(diag(1,3),30,3,2)
covu<-diag(1,3)
Sp=diag(0,3);
Sup=diag(1,3);
a=5
nburn=as.integer(10);
nbetween=as.integer(10);
nimp=as.integer(2);

# And finally we can run either the model with fixed or random cluster-specific covariance matrices:

imp<-jomo1rancathr(Y_cat, Y_numcat, X,Z,clus,betap,up,covp, covu,Sp,Sup,nburn,nbetween,nimp, a)
}
