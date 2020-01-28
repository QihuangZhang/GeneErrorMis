\name{GEEint}
\alias{GEEint}
\title{
GEE approach with internal validation data
}
\description{
The function implements the analysis of data with mixed type reponse using GEE approach. The continous response is subject to measurement error and the discrete function is subject to misclassification. In this scenario the mismeasurement is unknown, but we can estimate it by internal validation data.
}
\usage{
GEEint(initial, data.validation, data.mismeasure, 
                   dmindex1, dmindex2, vdindex1, vdindex2, 
                   CMval1index, CMval2index, CMmea1index, CMmea2index)
}
\arguments{
  \item{initial}{
a vector specifying the initial values of the regression parameters. Usually taking the value of the naive model. See example.
}
  \item{data.validation}{
a data.frame of the validation data. The validation data included the error-prone responses, true responses, covariates and the design matrix for the measurement error and misclassification models.
}
  \item{data.mismeasure}{
a data.frame of the non-validation data. The non-validation data included the error-prone responses, covariates and the design matrix for the measurement error and misclassification models. The precise measurement responses are not available.
}
  \item{dmindex1}{
a vector of the indices of the design matrix of the first response model in the non-validation data.
}
  \item{dmindex2}{
a vector of the indices of the design matrix of the second response model in the non-validation data.
}
  \item{vdindex1}{
a vector of the indices of the design matrix of the first response model in the validation data.
}
  \item{vdindex2}{
a vector of the indices of the design matrix of the second response model in the validation data.
}
  \item{CMval1index}{
a vector of the indices of the design matrix of the measurement error model in the validation data.
}
  \item{CMval2index}{
a vector of the indices of the design matrix of the misclassification model in the validation data.
}
  \item{CMmea1index}{
a vector of the indices of the design matrix of the measurement error model in the non-validation data.
}
  \item{CMmea2index}{
a vector of the indices of the design matrix of the misclassification model in the non-validation data.
}
}
\details{
The specification of \code{gamma} should correspond to the specification of the design matrix of \code{CovMis1}. The number of columns of \code{CovMis1} should be greater or equal to 2. The first column corresponds to the intercept. The user can take the vector to be all 1. The second column corresponds to the effect of the binary response. The users can specify it to be zero to indicate the second response does not affect the measurement error model.
}
\value{
\item{coefficients}{the coefficient of the main model after correction.}
\item{vcov}{a adjusted variance-covariance matrix estimated by sandwich estimator.}
\item{...}{Other components are the arguments that have been used in the function.}
}

\references{Grace Y Y. Statistical Analysis with Measurement Error or Misclassification[M]. Springer-Verlag New York, 2017.}

\author{
Qihuang Zhang and Grace Y. Yi
}


\seealso{
\code{\link[stats]{glm}}
}
\examples{
## load data
data(ToyGEEknown)

## Naive model
naive.model1 <- lm(Y1star ~ X + W, data = ToyGEEknown$naivedata)
naive.model2 <- glm(Y2star ~ X + W, family = binomial(link = logit), 
                    data = ToyGEEknown$naivedata)

intial <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)

## proposed method
example_known <- GEEknown(intial, Y1star = ToyGEEknown$naivedata$Y1star, 
                          Y2star = ToyGEEknown$naivedata$Y2star, 
                          DesignMatrix1 = ToyGEEknown$DesignMatrix, 
                          DesignMatrix2 = ToyGEEknown$DesignMatrix, 
                          CovMis1 = ToyGEEknown$CovMis1, 
                          CovMis2 = ToyGEEknown$CovMis2, 
                          gamma1=1, gamma=c(0,0.5), alpha1=-2, alpha0=-2, 
                          sigma_e=0.1)
example_known
}
\keyword{funtion}
\keyword{core}
\keyword{GEE}
