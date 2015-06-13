\name{jomo1con}
\alias{jomo1con}
 
\title{
JM Imputation of single level data with continuous variables only
}
\description{
Impute a single level dataset with continuous outcomes only. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler. Categorical covariates may be considered, but they have to be included with dummy variables.
}
\usage{
jomo1con(Y, X=matrix(1,nrow(Y),1), betap=matrix(0,ncol(X),ncol(Y)), 
covp=diag(1,ncol(Y)), Sp=diag(1,ncol(Y)), nburn=100, nbetween=100, nimp=5)
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
Number of burn in iterations. Default is 100.
}
  \item{nbetween}{
Number of iterations between two successive imputations. Default is 100.
}
  \item{nimp}{
Number of Imputations. Default is 5.
}

}
\details{
The Gibbs sampler algorithm used is described in detail in Chapter 3 of Carpenter and Kenward (2013). Regarding the choice of the priors, a flat prior is considered for beta, while an inverse-Wishart prior is given to the covariance matrix, with p-1 degrees of freedom, aka the minimum possible, to guarantee the greatest uncertainty. Binary or continuous covariates in the imputation model may be considered without any problem, but when considering a categorical covariate it has to be included through dummy variables (binary indicators) only. 
}
\value{
On screen, the posterior mean of the fixed effects estimates and of the covariance matrix are shown. The only argument returned is the imputed dataset in long format. Column "Imputation" indexes the imputations. Imputation number 0 are the original data.
}
\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapter 3, Wiley, ISBN: 978-0-470-74052-1.
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
nbetween=as.integer(200);
nimp=as.integer(5);

# Then we run he function:

imp<-jomo1con(Y,X,betap,covp,Sp,nburn,nbetween,nimp)

#Finally we can analyse the imputed datasets and combine results:

estimates<-rep(0,5)
ses<-rep(0,5)
estimates2<-rep(0,5)
ses2<-rep(0,5)
for (i in 1:5) {
  dat<-imp[imp$Imputation==i,]
  fit<-lm(measure~age+sex,data=dat)
  estimates[i]<-coef(summary(fit))[2,1]
  ses[i]<-coef(summary(fit))[2,2]
  estimates2[i]<-coef(summary(fit))[3,1]
  ses2[i]<-coef(summary(fit))[3,2]
}

#Here we use the BaBooN library to apply Rubins rules:

#library("BaBooN")
#MI.inference(estimates, ses^2)
#MI.inference(estimates2, ses2^2)

}
