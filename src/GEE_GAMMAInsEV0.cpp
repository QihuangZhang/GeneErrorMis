# include <Rcpp.h>
# include<iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix GEE_GAMMAInsEV0(NumericVector Y1star, NumericVector Y2star, NumericVector Y1, NumericVector Y2,
                             NumericMatrix CovMis1, NumericMatrix CovMis2,  int ncov1, int ncov2,
                                    double gamma1, NumericVector gamma, NumericVector alpha1, NumericVector alpha0, double sigma_e, 
                                    int fixgamma1, NumericVector fixgamma, int fixsigma_e, NumericVector fixalpha1, NumericVector fixalpha0){
  int i,j,k;
  int nval = Y1star.size();
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
  
  NumericMatrix GammaMatrix0(ncov1+ncov2+2+nq,ncov1+ncov2+2+nq);
  NumericVector mta_1c(nval),mta_20(nval),mta_21(nval);
  double pi0,pi1;
  int index;
  
  // Evaluate GammaMatrixI
  for (j = 0; j < ncov1 + ncov2 + 2 + nq; ++j){
    for (k = 0; k < ncov1 + ncov2 + 2 + nq; ++k){
      GammaMatrix0(j,k) = 0;
    }
  }
  
  for(i = 0; i < nval; ++i){ 
    for (j = 0; j < ncovMis; ++j){
      mta_21(i) += alpha1(j) * CovMis2(i,j);
      mta_20(i) += alpha0(j) * CovMis2(i,j);
    }
  }
  
  for(i = 0; i < nval; ++i){
    pi0 = exp(mta_20(i))/(1+exp(mta_20(i)));
    pi1 = exp(mta_21(i))/(1+exp(mta_21(i)));
    if (pi0 != pi0) {pi0 = 1;}
    if (pi1 != pi1) {pi1 = 1;}
    
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
