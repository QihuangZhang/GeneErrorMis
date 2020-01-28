# include <Rcpp.h>
# include<iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix GEE_GAMMAInsIV0(NumericVector Y1star, NumericVector Y2star, NumericVector Y1, NumericVector Y2,
                             NumericMatrix CovMis1, NumericMatrix CovMis2, NumericMatrix ValidationMatrix1, NumericMatrix ValidationMatrix2,
                                    NumericVector beta1, NumericVector beta2, double xi, double sigma,
                                    double gamma1, NumericVector gamma, NumericVector alpha1, NumericVector alpha0, double sigma_e, 
                                    int fixgamma1, NumericVector fixgamma, int fixsigma_e, NumericVector fixalpha1, NumericVector fixalpha0){
  int i,j,k,t;
  int nval = ValidationMatrix1.nrow();
  int ncov1 = ValidationMatrix1.ncol();
  int ncov2 = ValidationMatrix2.ncol();
  int ncovMea = CovMis1.ncol();
  int ncovMis = CovMis2.ncol();
  
  
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
  NumericMatrix Score(nval,nq);
  
  NumericMatrix Vii(2,2),Viinverse(2,2);
  NumericMatrix Dmatrix1(ncov1+ncov2+2,ncov1+ncov2+2),Dmatrix2(ncov1+ncov2+2,ncov1+ncov2+2);
  NumericMatrix GammaMatrix0(ncov1+ncov2+2+nq,ncov1+ncov2+2+nq);
  NumericVector eta1(nval),eta2(nval),eta2expit(nval);
  NumericVector mta_1c(nval),mta_20(nval),mta_21(nval);
  NumericVector vi2(nval);
  double pi0,pi1;
  int index;
  
  sigma = fabs(sigma);
  
  // Evaluate GammaMatrixI
  for (j = 0; j < ncov1 + ncov2 + 2 + nq; ++j){
    for (k = 0; k < ncov1 + ncov2 + 2 + nq; ++k){
      GammaMatrix0(j,k) = 0;
    }
  }
  
  for (j = 0; j < ncov1 + ncov2 + 2; ++j){
    for (k = 0; k < ncov1 + ncov2 + 2; ++k){
      Dmatrix1(j,k) = 0;
      Dmatrix2(j,k) = 0;
    }
  }
  
  for(i = 0; i < nval; ++i){ 
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
  }
  
  for(i = 0; i < nval; ++i){
    vi2(i) = eta2expit(i) * (1 - eta2expit(i));
    
    Vii(0,0) = pow(sigma,2);
    Vii(0,1) = xi*sigma*pow(vi2(i),0.5);
    Vii(1,0) = Vii(0,1);
    Vii(1,1) = vi2(i);
    
    Viinverse(0,0) = 1 / pow(sigma,2) / ( 1 - pow(xi,2)); 
    Viinverse(0,1) = - xi / (1 - pow(xi,2)) / sigma / pow(vi2(i),0.5); 
    Viinverse(1,0) = Viinverse(0,1);
    Viinverse(1,1) = 1 / (1 - pow(xi,2)) / vi2(i);
    
    for(j = 0; j < ncov1; ++j){
      Dmatrix1(j,0) = ValidationMatrix1(i,j) * Viinverse(0,0);
      Dmatrix1(j,1) = ValidationMatrix1(i,j) * Viinverse(0,1);
    }
    for(j = ncov1; j < ncov1+ncov2; ++j){
      Dmatrix1(j,0) = ValidationMatrix2(i,j-ncov1) * vi2(i) * Viinverse(1,0);
      Dmatrix1(j,1) = ValidationMatrix2(i,j-ncov1) * vi2(i) * Viinverse(1,1);
    }
    
    pi0 = exp(mta_20(i))/(1+exp(mta_20(i)));
    pi1 = exp(mta_21(i))/(1+exp(mta_21(i)));
    if (pi0 != pi0) {pi0 = 1;}
    if (pi1 != pi1) {pi1 = 1;}
    
    Dmatrix1(ncov1+ncov2,2) = 2 * sigma;
    Dmatrix1(ncov1+ncov2,3) = xi * pow(vi2(i),0.5);
    Dmatrix1(ncov1+ncov2,4) = xi * pow(vi2(i),0.5);
    
    Dmatrix1(ncov1+ncov2+1,2) = 0;
    Dmatrix1(ncov1+ncov2+1,3) = sigma  * pow(vi2(i),0.5);
    Dmatrix1(ncov1+ncov2+1,4) = sigma  * pow(vi2(i),0.5);
    Dmatrix1(ncov1+ncov2+1,5) = 0;
    
    // Construct the Dmatrix2
    for(j = 0; j < ncov1; ++j){
      Dmatrix2(0,j) = - ValidationMatrix1(i,j);
    }
    
    for(j = ncov1; j < ncov1+ncov2; ++j){
      Dmatrix2(1,j) = - ValidationMatrix2(i,j-ncov1) * vi2(i);
    }
    
    for(j = ncov1; j < ncov1 + ncov2; ++j){
      Dmatrix2(3,j) = - 0.5 * xi * sigma * (1-2 * eta2expit(i)) * ValidationMatrix2(i,j-ncov1) * pow(vi2(i),0.5);
      Dmatrix2(4,j) = Dmatrix2(3,j);
      Dmatrix2(5,j) = - (1-2 * eta2expit(i)) * ValidationMatrix2(i,j-ncov1) * vi2(i);
    }
    
    Dmatrix2(2,ncov1+ncov2) = - 2 * sigma;
    Dmatrix2(3,ncov1+ncov2) = - xi * pow(vi2(i),0.5);
    Dmatrix2(3,ncov1+ncov2+1) = - sigma * pow(vi2(i),0.5);
    Dmatrix2(4,ncov1+ncov2) = - xi * pow(vi2(i),0.5);
    Dmatrix2(4,ncov1+ncov2+1) = - sigma * pow(vi2(i),0.5);
    
    
    
    for(j = 0; j < ncov1 + ncov2 + 2; ++j)
    {
      for(k = 0; k < ncov1 + ncov2 + 2; ++k)
      {
        for(t=0; t<  ncov1 + ncov2 + 2; ++t)
        {
          GammaMatrix0(j,k) += Dmatrix1(j,t) * Dmatrix2(t,k);
        }
      }
    }
    
    index = -1;

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
      Score(i,index) = (Y1star(i)  - mta_1c(i)) * Y1(i)  / pow(sigma_e,2);
      }

    for (j = 0; j < fixgamma.length();++j){
      if (fixgamma(j)!=1) {
        index = index +1;
        if (j == 1) {
          Score(i,index) = ((Y1star(i)  - mta_1c(i)) * CovMis1(i,j) * Y2(i))  / pow(sigma_e,2);
         } else {
           Score(i,index) = ((Y1star(i)  - mta_1c(i)) * CovMis1(i,j))  / pow(sigma_e,2);
         }
        }
    }

    if (fixsigma_e!=1) {
      index = index +1;
      Score(i,index) =  -1/sigma_e +  pow((Y1star(i) - mta_1c(i)),2)/pow(sigma_e,3);
    }

    for (j = 0; j < fixalpha1.length();++j){
      if (fixalpha1(j)!=1) {
        index = index +1;
        Score(i,index) =  Y2(i) * CovMis2(i,j) * ((1-Y2star(i)) - pi1);
      }
    }

    for (j = 0; j < fixalpha0.length();++j){
      if (fixalpha0(j)!=1) {
        index = index +1;
        Score(i,index) =  (1-Y2(i)) * CovMis2(i,j) * (Y2star(i) - pi0);
      }
    }


    for(j = 0; j < nq; ++j)
    {
      for(k = 0; k < nq; ++k)
      { GammaMatrix0(ncov1 + ncov2 + 2 + j, ncov1 + ncov2 + 2 + k) += - Score(i,j) * Score(i,k);}
    }
    
    
  }
  
  return(GammaMatrix0);
}
                                         
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
