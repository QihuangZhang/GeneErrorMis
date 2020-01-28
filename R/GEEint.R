GEEint <- function(initial, data.validation, data.mismeasure, 
                   dmindex1, dmindex2, vdindex1, vdindex2, 
                   CMval1index, CMval2index, CMmea1index, CMmea2index)
 {model.measure <- lm(Y1star ~ -1 + offset(Y1) + Y2,data = data.validation) 
  model.class1 <- glm((1-Y2star) ~ 1, data = data.validation[data.validation$Y2==1,],family = binomial(link="logit")) 
  model.class0 <- glm(Y2star ~ 1, data = data.validation[data.validation$Y2==0,],family = binomial(link="logit")) 
 
  gamma2 <- model.measure$coefficients
  sigma_e <- sigma(model.measure)
  alpha1 <- model.class1$coefficients
  alpha0 <- model.class0$coefficients
 
  
  NR <- nleqslv(initial, GEE_UI_IV, jacobian=T, control=list(maxit=2000),
          data.validation = data.validation, data.mismeasure = data.mismeasure, 
          dmindex1=dmindex1, dmindex2=dmindex2, vdindex1=vdindex1, vdindex2=vdindex2, 
          CMmea1index=CMmea1index, CMmea2index=CMmea2index,
          gamma1 = 1, gamma=c(0,gamma2), alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)
  
  betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
 
  if (!any(is.na(betahat))) {
    cov <- GEE_covIV (betahat, data.validation, data.mismeasure,
                      dmindex1, dmindex2, vdindex1, vdindex2, 
                      CMval1index, CMval2index, CMmea1index, CMmea2index,
                      gamma1=1, gamma = c(0,gamma2), alpha1, alpha0, sigma_e,
                      fixgamma1=1, fixgamma=c(1,0), fixsigma_e=0, fixalpha1=0, fixalpha0=0)
    } else {
      cov <- NA
    }
 
 return(list(coefficients = betahat,
             vcov = cov))
} 