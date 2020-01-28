#### known mismeasurement mechanism ####

GEE_UI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                   CovMis1, CovMis2,
                   gamma1, gamma, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  ntheta <- length(theta)
  return(GEE_UfuncIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:((ntheta-2)/2)], theta[((ntheta-2)/2+1):(ntheta-2)],
                      sigma = theta[ntheta-1], xi = theta[ntheta], 
                      gamma1, gamma, alpha1=alpha1, alpha0=alpha0,
                      sigma_e))
}


GEE_SIGMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                      gamma1, gamma, alpha1, alpha0, sigma_e){
  ntheta <- length(theta)
  return(GEE_SIGMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                      theta[1:((ntheta-2)/2)], theta[((ntheta-2)/2+1):(ntheta-2)], 
                      sigma = theta[ntheta-1], xi = theta[ntheta],
                      gamma1, gamma, alpha1, alpha0, sigma_e))
}

GEE_GAMMA <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  ntheta <- length(theta)
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)], 
                        sigma = theta[ntheta-1], xi = theta[ntheta])
  return(GAMMA)
}

GEE_GAMMA.inv <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2){
  ntheta <- length(theta)
  GAMMA <- GEE_GAMMAIns(Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                        beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)],
                        sigma = theta[ntheta-1], xi = theta[ntheta])
  GAMMA.inv <- solve(GAMMA)
  return(GAMMA.inv)
}

GEE_cov <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                    gamma1, gamma, alpha1, alpha0, sigma_e){
  ntheta <- length(theta)
  GAMMA.inv <- GEE_GAMMA.inv(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2)
  SIGMA <- GEE_SIGMA(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                     gamma1, gamma, alpha1, alpha0, sigma_e)
  covmatrix <- GAMMA.inv %*% SIGMA %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
}

### internal validation ####

GEE_UI_IV <- function(theta,  data.validation, data.mismeasure, 
                      dmindex1, dmindex2, vdindex1, vdindex2, 
                      CMmea1index, CMmea2index,
                      gamma1, gamma, alpha1, alpha0, sigma_e){
  # cat(theta, " \n")
  ntheta <- length(theta)
  return(GEE_UfuncInsIV(Y1star=data.mismeasure$Y1star,
                        Y2star=data.mismeasure$Y2star,
                        Y1 = data.validation$Y1,
                        Y2 = data.validation$Y2,
                        DesignMatrix1 = as.matrix(data.mismeasure[,dmindex1]), 
                        DesignMatrix2 = as.matrix(data.mismeasure[,dmindex2]), 
                        ValidationMatrix1 = as.matrix(data.validation[,vdindex1]),
                        ValidationMatrix2 = as.matrix(data.validation[,vdindex2]),
                        CovMis1 = as.matrix(data.mismeasure[,CMmea1index]),
                        CovMis2 = as.matrix(data.mismeasure[,CMmea2index]),
                        beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)], sigma = theta[ntheta-1], xi = theta[ntheta], 
                        gamma1=gamma1, gamma=gamma, alpha1=alpha1, alpha0=alpha0, sigma_e=sigma_e
  ))
}



GEE_GAMMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  ntheta <- length(theta)
  return(GEE_GAMMAInsIV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ValidationMatrix1, ValidationMatrix2,
                         beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)], xi=theta[ntheta], sigma=theta[ntheta-1],
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}


GEE_GAMMA_IVI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  ntheta <- length(theta)
  return(GEE_GAMMAInsIVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)], xi=theta[ntheta], sigma=theta[ntheta-1],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_SIGMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  ntheta <- length(theta)
  return(GEE_GAMMAInsIV0(Y1star, Y2star, Y1, Y2,
                         CovMis1, CovMis2, ValidationMatrix1, ValidationMatrix2,
                         beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)], xi=theta[ntheta], sigma=theta[ntheta-1],
                         gamma1, gamma, alpha1, alpha0, sigma_e, 
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_SIGMA_IV0 <- function(theta, Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  ntheta <- length(theta)
  return(GEE_SIGMAInsIV0(Y1star, Y2star, Y1, Y2, ValidationMatrix1, ValidationMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)], xi=theta[ntheta], sigma=theta[ntheta-1],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}


GEE_SIGMA_IVI <- function(theta, Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2, 
                          gamma1, gamma, alpha1, alpha0, sigma_e,
                          fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  ntheta <- length(theta)
  return(GEE_SIGMAInsIVI(Y1star, Y2star, DesignMatrix1, DesignMatrix2,
                         CovMis1, CovMis2,
                         beta1=theta[1:((ntheta-2)/2)], beta2=theta[((ntheta-2)/2+1):(ntheta-2)], xi=theta[ntheta], sigma=theta[ntheta-1],
                         gamma1, gamma, alpha1, alpha0, sigma_e,
                         fixgamma1=fixgamma1, fixgamma=fixgamma, fixsigma_e=fixsigma_e, fixalpha1=fixalpha1, fixalpha0=fixalpha0)
  )
}

GEE_covIV <- function(theta, data.validation, data.mismeasure,
                      dmindex1, dmindex2, vdindex1, vdindex2, 
                      CMval1index, CMval2index, CMmea1index, CMmea2index,
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0){
  nvalidation <- dim(data.validation)[1]
  nsample <- dim(data.mismeasure)[1] + nvalidation
  
  M0 <- GEE_GAMMA_IV0(theta, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      ValidationMatrix1 = as.matrix(data.validation[,vdindex1]),
                      ValidationMatrix2 = as.matrix(data.validation[,vdindex2]), 
                      CovMis1 = as.matrix(data.validation[,CMval1index]), 
                      CovMis2 = as.matrix(data.validation[,CMval2index]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  M1 <- GEE_GAMMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,dmindex1]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,dmindex2]), 
                      CovMis1 = as.matrix(data.mismeasure[,CMmea1index]), 
                      CovMis2 = as.matrix(data.mismeasure[,CMmea2index]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  GAMMA_IV <- M1 + M0
  
  B0 <- GEE_SIGMA_IV0(theta, 
                      Y1star=data.validation$Y1star, 
                      Y2star=data.validation$Y2star, 
                      Y1 = data.validation$Y1, 
                      Y2 = data.validation$Y2, 
                      ValidationMatrix1 = as.matrix(data.validation[,vdindex1]),
                      ValidationMatrix2 = as.matrix(data.validation[,vdindex2]), 
                      CovMis1 = as.matrix(data.validation[,CMval1index]), 
                      CovMis2 = as.matrix(data.validation[,CMval2index]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  B1 <- GEE_SIGMA_IVI(theta, 
                      Y1star=data.mismeasure$Y1star, 
                      Y2star=data.mismeasure$Y2star, 
                      DesignMatrix1 = as.matrix(data.mismeasure[,dmindex1]),
                      DesignMatrix2 = as.matrix(data.mismeasure[,dmindex2]), 
                      CovMis1 = as.matrix(data.mismeasure[,CMmea1index]), 
                      CovMis2 = as.matrix(data.mismeasure[,CMmea2index]), 
                      gamma1, gamma, alpha1, alpha0, sigma_e,
                      fixgamma1,fixgamma,fixsigma_e,fixalpha1,fixalpha0)
  
  SIGMA_IV <- B1 + B0
  
  GAMMA.inv <- solve(GAMMA_IV,tol=1e-200)
  covmatrix <- GAMMA.inv %*% SIGMA_IV %*% t(as.matrix(GAMMA.inv))
  return(covmatrix)
  # return(list(M1=M1,M0=M0,B1=B1,B0=B0,GAMMA_IV=GAMMA_IV,SIGMA_IV=SIGMA_IV,covmatrix))
}