\name{cldata}
\alias{cldata}
\docType{data}
\title{
A simulated clustered dataset 
}
\description{
A simulated dataset to test functions for imputation of clustered data. 
}
\usage{data(cldata)}
\format{
  A data frame with 1000 observations on the following 6 variables.
  \describe{
    \item{\code{age}}{A numeric variable with individual ages. Fully observed.}
    \item{\code{measure}}{A numeric variable with some measure of interest (unspecified). This is partially observed.}
    \item{\code{sex}}{A binary variable with sex indicator. Fully observed.}
    \item{\code{social}}{A 4-category variable with some social status indicator. This is partially observed.}
    \item{\code{city}}{The cluster indicator vector. 10 cities are indexed 0 to 9.}
    \item{\code{id}}{The id for individuals within each city.}
  }
}
\details{
These are not real data, they are simulated to illustrate the use of the main functions of the package.}


