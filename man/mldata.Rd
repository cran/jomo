\name{mldata}
\alias{mldata}
\docType{data}
\title{
A simulated multilevel dataset 
}
\description{
A simulated dataset to test multilevel functions, i.e. jomo1rancon, jomo1rancat, jomo1ranmix, jomo1ranconhr,jomo1rancathr and jomoran1mixhr. 
}
\usage{data(mldata)}
\format{
  A data frame with 1000 observations on the following 6 variables.
  \describe{
    \item{\code{age}}{a numeric vector of ages. Fully observed.}
    \item{\code{measure}}{a numeric vector for our measure of interest (unspecified). This is partially observed.}
    \item{\code{sex}}{a binary vector for sex. Fully observed.}
    \item{\code{social}}{a 4-category vector for social status indicators. This is partially observed.}
    \item{\code{city}}{The cluster indicator vector. 10 cities are indexed 0 to 9.}
    \item{\code{id}}{The id within each city.}
  }
}
\details{
These are not real data, they are simulated to illustrate the use of the main functions of the package.}


