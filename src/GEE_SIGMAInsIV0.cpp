# include <Rcpp.h>
# include<iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix GEE_SIGMAInsIV0(NumericVector Y1star, NumericVector Y2star, NumericVector Y1, NumericVector Y2,
                              NumericMatrix ValidationMatrix1, NumericMatrix ValidationMatrix2, NumericMatrix CovMis1, NumericMatrix CovMis2, 
                              NumericVector beta1, NumericVector beta2, double xi, 
                              double sigma, double gamma1, NumericVector gamma, NumericVector alpha1, NumericVector alpha0, double sigma_e,
                              int fixgamma1, NumericVector fixgamma, int fixsigma_e, NumericVector fixalpha1, NumericVector fixalpha0){
  int i,j,k;
  int nsam = ValidationMatrix1.nrow();
  int ncov1 = ValidationMatrix1.ncol();
  int ncov2 = ValidationMatrix2.ncol();
  int ncovMea = CovMis1.ncol();
  int ncovMis = CovMis2.ncol();
  NumericMatrix Vii(2,2),Viinverse(2,2);
  NumericVector eta1(nsam),eta2(nsam),eta2expit(nsam);
  NumericVector mta_20(nsam),mta_21(nsam),mta_1c(nsam);
  NumericVector vi2(nsam);
  double pi0,pi1;
  double S1,S2,S3;

 
  sigma_e = fabs(sigma_e);
  sigma = fabs(sigma);
  
  
  // Dimension of the mismeasurement models
  int nq=1-fixsigma_e;
  for (i = 0; i < fixalpha1.length();++i){
    if (fixalpha1(i)!=1) {nq = nq +1;}
  }
  for (i = 0; i < fixalpha0.length();++i){
    if (fixalpha0(i)!=1) {nq = nq +1;}
  } 
  if (fixgamma1!=1) {nq = nq +1;}
  for (i = 0; i < fixgamma.length();++i){
    if (fixgamma(i)!=1) {nq = nq +1;}
  }
  
  NumericMatrix SIGMA(ncov1+ncov2+2+nq,ncov1+ncov2+2+nq);
  NumericMatrix Ufunc(nsam,ncov1+ncov2+2+nq);
  
  
  for(i = 0; i < nsam; ++i){ 
    eta1(i) = 0;
    eta2(i) = 0;
    for (j = 0; j < ncov1; ++j){
      eta1(i) += beta1(j) * ValidationMatrix1(i,j);
    }
    for (j = 0; j < ncov2; ++j){
      eta2(i) += beta2(j) * ValidationMatrix2(i,j);
    }
    eta2expit(i) = exp(eta2(i))/(1+exp(eta2(i)));
    if (eta2expit(i) != eta2expit(i)) {eta2expit(i) = 1;}
    for (j = 0; j < ncovMis; ++j){
      mta_21(i) += alpha1(j) * CovMis2(i,j);
      mta_20(i) += alpha0(j) * CovMis2(i,j);
    }
    for (j = 0; j < ncov1 + ncov2 + 2; ++j){
      Ufunc(i,j) = 0;
    }
  }
  
  for(i = 0; i < ncov1+ncov2+2; ++i){ 
    for(j = 0; j < ncov1+ncov2+2; ++j){ 
      SIGMA(i,j) = 0;
    }
  }
  
  for(i = 0; i < nsam; ++i){
    vi2(i) = eta2expit(i) * (1 - eta2expit(i));
    
    Vii(0,0) = pow(sigma,2);
    Vii(0,1) = xi*sigma*pow(vi2(i),0.5);
    Vii(1,0) = Vii(0,1);
    Vii(1,1) = vi2(i);
    
    Viinverse(0,0) = 1 / pow(sigma,2) / ( 1 - pow(xi,2)); 
    Viinverse(0,1) = - xi / (1 - pow(xi,2)) / sigma / pow(vi2(i),0.5); 
    Viinverse(1,0) = Viinverse(0,1);
    Viinverse(1,1) = 1 / (1 - pow(xi,2)) / vi2(i);
    
    
    pi0 = exp(mta_20(i))/(1+exp(mta_20(i)));
    pi1 = exp(mta_21(i))/(1+exp(mta_21(i)));
    if (pi0 != pi0) {pi0 = 1;}
    if (pi1 != pi1) {pi1 = 1;}
    
    S1 = pow(Y1(i)-eta1(i),2);
    S2 = (Y1(i)-eta1(i))*(Y2(i)-eta2expit(i)); 
    S3 = pow(Y2(i)-eta2expit(i),2);
    
    for (j = 0; j < ncov1; ++j){
      Ufunc(i,j) += ValidationMatrix1(i,j) * Viinverse(0,0) * (Y1(i) - eta1(i))
                  + ValidationMatrix1(i,j) * Viinverse(0,1) * (Y2(i) - eta2expit(i));
    }
    
    for (j = ncov1; j < ncov1 + ncov2; ++j){
      Ufunc(i,j) += ValidationMatrix2(i,j-ncov1) * eta2expit(i) * (1-eta2expit(i)) * Viinverse(1,0) * (Y1(i) - eta1(i))
      + ValidationMatrix2(i,j-ncov1) * eta2expit(i) * (1-eta2expit(i)) * Viinverse(1,1) * (Y2(i) - eta2expit(i));
    }
    
    Ufunc(i,ncov1 + ncov2) += 2 * sigma * (S1-Vii(0,0)) 
      + 2*xi * pow(vi2(i),0.5) * (S2-Vii(0,1));

    Ufunc(i,ncov1 + ncov2 + 1) += 2 * sigma * pow(vi2(i),0.5) * (S2-Vii(0,1));
                   
    int index = ncov1 + ncov2 + 1;
                   
    for (j = 0; j < ncovMea; ++j){
      if (j == 1) {
         mta_1c(i) += gamma(j) * CovMis1(i,j) * Y2(i);
      } else {
         mta_1c(i) += gamma(j) * CovMis1(i,j);
      }
    }
                   
    mta_1c(i) += gamma1 * Y1(i);
                   
    if (fixgamma1!=1) {
       index = index +1;
       Ufunc(i,index) = (Y1star(i)  - mta_1c(i)) * Y1(i)  / pow(sigma_e,2);
    }
                   
    for (j = 0; j < fixgamma.length();++j){
       if (fixgamma(j)!=1) {
         index = index +1;
            if (j == 1) {
              Ufunc(i,index) = ((Y1star(i)  - mta_1c(i)) * CovMis1(i,j) * Y2(i))  / pow(sigma_e,2);
            } else {
              Ufunc(i,index) = ((Y1star(i)  - mta_1c(i)) * CovMis1(i,j))  / pow(sigma_e,2);
              }
        }
    }
                   
    if (fixsigma_e!=1) {
       index = index +1;
       Ufunc(i,index) =  -1/sigma_e +  pow((Y1star(i) - mta_1c(i)),2)/pow(sigma_e,3);
    }
                   
    for (j = 0; j < fixalpha1.length();++j){
       if (fixalpha1(j)!=1) {
          index = index +1;
          Ufunc(i,index) =  Y2(i) * CovMis2(i,j) * ((1-Y2star(i)) - pi1);
       }
    }
                   
    for (j = 0; j < fixalpha0.length();++j){
       if (fixalpha0(j)!=1) {
          index = index +1;
          Ufunc(i,index) =  (1-Y2(i)) * CovMis2(i,j) * (Y2star(i) - pi0);
       }
    }                   
                   
    for(j = 0; j < ncov1+ncov2+2+nq; ++j){ 
       for(k = 0; k < ncov1+ncov2+2+nq; ++k){ 
         SIGMA(j,k) += Ufunc(i,j) * Ufunc(i,k);
      }
   }  
  }
  
  return(SIGMA);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
