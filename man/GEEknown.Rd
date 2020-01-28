\name{GEEknown}
\alias{GEEknown}
\title{
GEE approach with Known Mismeasurement
}
\description{
The function implements the analysis of data with mixed type reponse using GEE approach. The continous response is subject to measurement error and the discrete function is subject to misclassification. In this scenario the mismeasurement is known.
}
\usage{
GEEknown(initial, Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
         CovMis1, CovMis2, gamma1=1, gamma, alpha1, alpha0, sigma_e)
}
\arguments{
  \item{initial}{
a vector specifying the initial values of the regression parameters. Usually taking the value of the naive model. See example.
}
  \item{Y1star}{
a vector specifying the error-prone continouse response.
}
  \item{Y2star}{
a vector specifying the binary response subject to misclassification.
}
  \item{DesignMatrix1}{
a matrix of the design matrix indecating the covariates for the continous model.
}
  \item{DesignMatrix2}{
a matrix of the design matrix indecating the covariates for the binary model.
}
  \item{CovMis1}{
a matrix of the design matrix for the measurement error model. The first two column are fixed for the intercept and the effet the second response. 
}
  \item{CovMis2}{
a matrix of the design matrix for the misclassification models. We assume the structure are the same for both misclassification models.
}
  \item{gamma1}{
a scaler indecating the value of the effect of the true continouse response in the measurement error model. The initial value is taken as 1.
}
  \item{gamma}{
a numeric vector specifying the values of other parameters for the measurement error model. 
}
  \item{alpha1}{
a vector of the parameters of misclassification model when the true binary response equals to 1.
}
  \item{alpha0}{
a vector of the parameters of misclassification model when the true binary response equals to 0.
}
  \item{sigma_e}{
a scaler indecating the standard error parameters for the measurement error model.
}
}
\details{
The specification of \code{gamma} should corresponds to the specification of the design matrix of \code{CovMis1}. The number of columns of \code{CovMis1} should be greater or equal to 2. The first column corresponds to the intercept. The user can take the vector to be all 1. The second column corresponds to the effect of the binary response. The users can sepecify it to be zero to indecate the second response does not effect the measurement error model.
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

initial <- c(naive.model1$coefficients,naive.model2$coefficients,1,0)

## proposed method
example_known <- GEEknown(initial, Y1star = ToyGEEknown$naivedata$Y1star, 
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
