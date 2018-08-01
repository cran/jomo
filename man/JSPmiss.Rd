\name{JSPmiss}
\alias{JSPmiss}
\docType{data}
\title{
Exam results for six inner London Education Authorities
}
\description{
A partially observed version of the jspmix1 dataset in package R2MLwiN. This is an educational dataset of pupils' test scores, a subset of the Junior School Project (Mortimore et al, 1988).
}
\usage{data(cldata)}
\format{
  A data frame with 4059 observations on the following 6 variables.
  \describe{
    \item{\code{school}}{A school identifier.}
    \item{\code{id}}{A student ID.}
    \item{\code{fluent}}{Fluency in English indicator, where 0 = beginner, 1 = intermediate, 2 = fully fluent; measured in Year 1.}
    \item{\code{sex}}{Sex of pupil; numeric with levels 0 (boy), 1 (girl).}
    \item{\code{cons}}{A column of 1s. Useful to add an intercept to th eimputation model.}
    \item{\code{ravens}}{Test score, out of 40; measured in Year 1.}
    \item{\code{english}}{Pupils' English test score, out of 100; measured in Year 3.}
    \item{\code{behaviour}}{Pupils' behaviour score, where lowerquarter = pupil rated in bottom 25\%, and upper otherwise; measured in Year 3.}

  }
}
\details{
These fully observed verison of the data is available with package R2MLwiN.}
\source{
Browne, W. J. (2012) MCMC Estimation in MLwiN Version 2.26. University of Bristol: Centre for Multilevel Modelling.

Mortimore, P., Sammons, P., Stoll, L., Lewis, D., Ecob, R. (1988) School Matters. Wells: Open Books.

Rasbash, J., Charlton, C., Browne, W.J., Healy, M. and Cameron, B. (2009) MLwiN Version 2.1. Centre for Multilevel Modelling, University of Bristol.


}

