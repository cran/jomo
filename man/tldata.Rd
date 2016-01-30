\name{tldata}
\alias{tldata}
\docType{data}
\title{
A simulated 2-level dataset 
}
\description{
A simulated dataset to test 2-level functions, i.e. jomo2com and jomo2hr. 
}
\usage{data(tldata)}
\format{
  A data frame with 1000 observations on the following 6 variables.
  \describe{
    \item{\code{measure.a}}{A numeric variable with some measure of interest (unspecified). This is partially observed.}
    \item{\code{measure.b}}{A numeric variable with some measure of interest (unspecified). This is fully observed.}
      \item{\code{measure.a2}}{A numeric variable with some level-2 measure of interest (unspecified). This is partially observed.}
    \item{\code{previous.events}}{A binary variable indicating if a patient has previous history of (unspecified) events. Patially observed.}
    \item{\code{group}}{A 3-category variable indicating to which group each patient belongs. This is partially observed.}
     \item{\code{big.city}}{A binary variable indicating if each city has more than 100000 inhabitants. Patially observed.}
    \item{\code{region}}{A 3-category variable indicating to which region each city belongs. This is fully observed.}
    \item{\code{city}}{The cluster indicator vector. 200 cities are indexed 0 to 199.}
    \item{\code{id}}{The id for each individual within each city.}
  }
}
\details{
These are not real data, they are simulated to illustrate the use of the main functions of the package.}


