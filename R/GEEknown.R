GEEknown <- function(initial, Y1star, Y2star, DesignMatrix1, DesignMatrix2, 
                     CovMis1, CovMis2, gamma1 = 1, gamma, alpha1, alpha0, sigma_e)
 { 
  NR <- nleqslv(x=initial, fn=GEE_UI, Y1star=Y1star, Y2star=Y2star, DesignMatrix1=DesignMatrix1, DesignMatrix2=DesignMatrix2,
                 CovMis1=CovMis1, CovMis2=CovMis2, control=list(maxit=2000),
                 gamma1 = gamma1, gamma=gamma, alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)
 
 betahat <- ifelse(abs(NR$x)<10,NR$x,NA)
 
 if (!any(is.na(betahat))) {
   cov <- GEE_cov(betahat,Y1star, Y2star, DesignMatrix1, DesignMatrix2, CovMis1, CovMis2,
                  gamma1 = gamma1, gamma=gamma, alpha1= alpha1, alpha0= alpha0, sigma_e = sigma_e)
   } else {
     cov <- NULL
   }
 
 return(list(coefficients = betahat,
             vcov = cov))
} 