\name{surdata}
\alias{surdata}
\docType{data}
\title{
  A simulated dataset with survival data
}
\description{
  A simulated dataset to test functions for imputation compatible with cox model.
}
\usage{data(cldata)}
\format{
  A data frame with 500 observations on the following 5 variables.
  \describe{
    \item{\code{measure}}{A numeric variable with some measure of interest (unspecified). This is partially observed.}
    \item{\code{sex}}{A binary variable with sex indicator. Partially observed.}
    \item{\code{id}}{The id for individuals within each city.}
    \item{\code{time}}{Time to event (death or censoring).}
    \item{\code{status}}{Binary variables, which takes value 0 for censored observations and 1 for deaths/events.}
  }
}
\details{
  These are not real data, they are simulated to illustrate the use of the main functions of the package.}


