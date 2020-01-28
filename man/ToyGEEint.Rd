\name{ToyGEEint}
\docType{data}
\alias{ToyGEEint}
\title{Toy example data with internal validation data}
\description{
  This data set gives an example data for the illustration of usage of \code{\link{GEEint}} function. The dataset contains validation data and non-validation. The non-validation data contains error-prone responses. Validation data have both error-prone and true responses. Both validation data and non-validation data contain the design matrix for the covariates in the response model, measurement error model and misclassification model.
}
\usage{data(ToyGEEint)}
\format{A list of two data frames. Non-validation data, a data frame 500 observations of 8 variables. Validation data, 500 observations of 10 variables. }
\keyword{datasets}