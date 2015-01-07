\name{jomo1ranmix}
\alias{jomo1ranmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JM Imputation of clustered data with mixed variable types
}
\description{
Impute a clustered dataset with mixed data types as outcome. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler where the covariance matrix is updated either with a Metropolis-Hastings step or from the inverse Wishart distribution. Fully observed categorical covariates may be considered as covariates as well, but they have to be included as dummy variables.
}
\usage{
jomo1ranmix(Y_con, Y_cat, Y_numcat, X=matrix(1,nrow(Y_cat),1), Z=matrix(1,nrow(Y_cat),1),
clus, betap=matrix(0,ncol(X),(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))), 
up=matrix(0,nrow(unique(clus)),ncol(Z)*(ncol(Y_con)+(sum(Y_numcat)-length(Y_numcat)))),
covp=diag(1,ncol(betap)), covu=diag(1,ncol(up)), Sp=diag(1,ncol(covp)), 
Sup=diag(1,ncol(covu)), nburn=100, nbetween=100, nimp=5, meth="MH") 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_con}{
A data frame, or matrix, with continuous responses of the joint imputation model. Rows correspond to different observations, while columns are different variables. Missing values are coded as NA. If no continuous outcomes are present in the model, the default for Y_con is NULL.
}
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
 \item{meth}{
Estimation method for covariance matrix. When meth="MH", the covariance matrix at each iteration is estimated with a Metropolis Hastings step, while when meth="IW" is specified, then covariance matrix is drawn from inverse Wishart distribution and afterwards constraints on elements related to latent normals are imposed.
}

}
\details{
The Gibbs sampler algorithm used is described in detail in Chapter 9 of Carpenter and Kenward (2013). Regarding the choice of the priors, a flat prior is considered for beta, while an inverse-Wishart prior is given to the covariance matrix, with p-1 degrees of freedom, aka the minimum possible, to guarantee the greatest uncertainty. A Metropolis Hastings step is implemented to update the covariance matrix when option meth="MH", as described in the book. When option meth="IW", then we update the covariance matrix drawing from the inverse Wishart distribution as we usually do when only continuous data are present in the model, and afterwards we constrain variances of the latent normals to 1 and covariances among latent normals related to the same categorical variable to 0. Then positive definitiveness of the matrix is checked and, if not met, another matrix is drawn.  Binary or continuous covariates in the imputation model may be considered without any problem, but when considering a categorical covariate it has to be included with dummy variables (binary indicators) only. 
}
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 9, Wiley, ISBN: 978-0-470-74052-1.
}

\examples{

#First of all we load and attach the data:

data(mldata)
attach(mldata)

#Then we define the inputs:
# nimp, nburn and nbetween are smaller than they should. This is
#just because of CRAN policies on the examples.

Y_con=data.frame(measure,age)
Y_cat=data.frame(social)
Y_numcat=matrix(4,1,1)
X=data.frame(rep(1,1000),sex)
Z<-data.frame(rep(1,1000))
clus<-data.frame(city)
betap<-matrix(0,2,5)
up<-matrix(0,10,5)
covp<-diag(1,5)
covu<-diag(1,5)
Sp=diag(0,5);
Sup=diag(1,5);
nburn=as.integer(100);
nbetween=as.integer(100);
nimp=as.integer(5);

#Then we can run the function, deciding which method we want to use to update the covariance matrix:

imp<-jomo1ranmix(Y_con, Y_cat, Y_numcat, X,Z,clus,betap,up,covp, covu,Sp,Sup,nburn,nbetween,nimp, 
                meth="MH")

# Or:
# imp<-jomo1ranmix(Y_con, Y_cat, Y_numcat, X,Z,clus,betap,up,covp, covu,Sp,Sup,nburn,nbetween,nimp, 
#                meth="IW")

# We run our substantive model on the 5 imputed datasets:

estimates<-matrix(0,5,5)
ses<-matrix(0,5,5)
for (i in 1:5) {
  dat<-imp[imp$Imputation==i,]
  fit<-lm(measure~age+sex+factor(social)+factor(clus),data=dat)
  estimates[i,1:5]<-coef(summary(fit))[2:6,1]
  ses[i,1:5]<-coef(summary(fit))[2:6,2]
}

# And finally we aggregate results with Rubin's rules, using BaBooN package:

#library("BaBooN")
#MI.inference(estimates[,1], ses[,1])
#MI.inference(estimates[,2], ses[,2])
#MI.inference(estimates[,3], ses[,3])
#MI.inference(estimates[,4], ses[,4])
#MI.inference(estimates[,5], ses[,5])
}
