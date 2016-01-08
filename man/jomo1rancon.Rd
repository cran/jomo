\name{jomo1rancon}
\alias{jomo1rancon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JM Imputation of clustered data with continuous variables only
}
\description{
Impute a clustered dataset with continuous outcomes only. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler. Categorical covariates may be considered, but they have to be included with dummy variables.
}
\usage{
jomo1rancon(Y, X=matrix(1,nrow(Y),1), Z=matrix(1,nrow(Y),1), clus, 
            beta.start=matrix(0,ncol(X),ncol(Y)),u.start=NULL, 
            l1cov.start=diag(1,ncol(Y)),l2cov.start=NULL, 
            l1cov.prior=diag(1,ncol(Y)), l2cov.prior=NULL,
            nburn=100, nbetween=100, nimp=5, output=1, out.iter=10)
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
  \item{beta.start}{
Starting value for beta, the vector(s) of fixed effects. Rows index different covariates and columns index different outcomes. The default is a matrix of zeros.
}
  \item{u.start}{
A matrix where different rows are the starting values within each cluster for the random effects estimates u. The default is a matrix of zeros.
}
  \item{l1cov.start}{
Starting value for the covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model. The default is the identity matrix.
}
  \item{l2cov.start}{
Starting value for the level 2 covariance matrix. Dimension of this square matrix is equal to the number of outcomes in the imputation model times the number of random effects. The default is an identity matrix.
}
  \item{l1cov.prior}{
Scale matrix for the inverse-Wishart prior for the covariance matrix. The default is the identity matrix.
}
  \item{l2cov.prior}{
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

 \item{output}{
When set to any value different from 1 (default), no output is shown on screen at the end of the process.
}
\item{out.iter}{
When set to K, every K iterations a message "Iteration number N*K completed" is printed on screen. Default is 10. 
}

}
\details{
The Gibbs sampler algorithm used is a simplification of the one described in detail in Chapter 9 of Carpenter and Kenward (2013), where we exclude the presence of level 2 variables. Regarding the choice of the priors, a flat prior is considered for beta, while an inverse-Wishart prior is given to the covariance matrices, with p-1 degrees of freedom, aka the minimum possible, to guarantee the greatest uncertainty. Binary or continuous covariates in the imputation model may be considered without any problem, but when considering a categorical covariate it has to be included with dummy variables (binary indicators) only. 
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

#Then we define all the inputs:
Y<-data.frame(measure,age)
X<-data.frame(rep(1,1000),sex)
Z<-data.frame(rep(1,1000))
clus<-data.frame(city)
beta.start<-matrix(0,2,2)
u.start<-matrix(0,10,2)
l1cov.start<-diag(1,2)
l2cov.start<-diag(1,2)
l1cov.prior=diag(1,2);
nburn=as.integer(200);
nbetween=as.integer(200);
nimp=as.integer(5);
l2cov.prior=diag(1,5);

#And finally we run the imputation function:
imp<-jomo1rancon(Y,X,Z,clus,beta.start,u.start,l1cov.start, l2cov.start,l1cov.prior,
             l2cov.prior,nburn,nbetween,nimp)

cat("Original value was missing(",imp[4,1],"), imputed value:", imp[1004,1])

#Then we run the model on the imputed datsets:

estimates<-rep(0,5)
ses<-rep(0,5)
estimates2<-rep(0,5)
ses2<-rep(0,5)
for (i in 1:5) {
  dat<-imp[imp$Imputation==i,]
  fit<-lm(measure~age+sex+factor(clus),data=dat)
  estimates[i]<-coef(summary(fit))[2,1]
  ses[i]<-coef(summary(fit))[2,2]
  estimates2[i]<-coef(summary(fit))[3,1]
  ses2[i]<-coef(summary(fit))[3,2]
}

#And finally we aggregate results through RUbin's rules with package BaBooN.

#library("BaBooN")
#MI.inference(estimates, ses^2)
#MI.inference(estimates2, ses2^2)

}
