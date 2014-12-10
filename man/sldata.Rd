\name{sldata}
\alias{sldata}
\docType{data}
\title{
A simulated single level dataset 
}
\description{
A simulated dataset to test single level functions, i.e. jomo1con, jomo1cat and jomo1mix. 
}
\usage{data(sldata)}
\format{
  A data frame with 300 observations on the following 4 variables.
  \describe{
    \item{\code{age}}{a numeric vector of ages. Fully observed.}
    \item{\code{measure}}{a numeric vector for our measure of interest (unspecified). This is partially observed.}
    \item{\code{sex}}{a binary vector for sex. Fully observed.}
    \item{\code{social}}{a 4-category vector for social status indicators. This is partially observed.}

  }
}
\details{
These are not real data, they are simulated to illustrate the use of the main functions of the package.}

