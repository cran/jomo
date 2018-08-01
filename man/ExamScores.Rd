\name{ExamScores}
\alias{ExamScores}
\docType{data}
\title{
Exam results for six inner London Education Authorities
}
\description{
A partially observed version of the tutorial dataset in package R2MLwiN.It includes examination results from six inner London Education Authorities (school boards).
}
\usage{data(cldata)}
\format{
  A data frame with 4059 observations on the following 6 variables.
  \describe{
    \item{\code{school}}{A school identifier.}
    \item{\code{student}}{A student ID.}
    \item{\code{normexam}}{Students' exam score at age 16, normalised and partially observed.}
    \item{\code{sex}}{Sex of pupil; a factor with levels boy, girl.}
    \item{\code{cons}}{A column of 1s. Useful to add an intercept to th eimputation model.}
    \item{\code{standlrt}}{Students' score at age 11 on the London Reading Test (LRT), standardised.}
    \item{\code{schgend}}{Schools' gender; a factor with levels corresponding to mixed school (mixedsch), boys' school (boysch), and girls' school (girlsch).}
    \item{\code{avslrt}}{Average LRT score in school.}
    \item{\code{schav}}{Average LRT score in school, coded into 3 categories: low = bottom 25\%, mid = middle 50\%, high = top 25\%.}
    \item{\code{vrband}}{Students' score in test of verbal reasoning at age 11, a factor with 3 levels: vb1 = top 25\%, vb2 = middle 50\%, vb3 = bottom 25\%.}

  }
}
\details{
These fully observed verison of the data is available with package R2MLwiN.}


\source{
Browne, W. J. (2012) MCMC Estimation in MLwiN Version 2.26. University of Bristol: Centre for Multilevel Modelling.

Goldstein, H., Rasbash, J., Yang, M., Woodhouse, G., Pan, H., Nuttall, D., Thomas, S. (1993) A multilevel analysis of school examination results. Oxford Review of Education, 19, 425-433.

Rasbash, J., Charlton, C., Browne, W.J., Healy, M. and Cameron, B. (2009) MLwiN Version 2.1. Centre for Multilevel Modelling, University of Bristol.


}