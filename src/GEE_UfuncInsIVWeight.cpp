# include <Rcpp.h>
# include<iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector GEE_UfuncInsIVWeight(NumericVector Y1star, NumericVector Y2star, NumericVector Y1, NumericVector Y2,
                            NumericMatrix DesignMatrix1, NumericMatrix DesignMatrix2, NumericMatrix ValidationMatrix1, NumericMatrix ValidationMatrix2,
                            NumericMatrix CovMis1, NumericMatrix CovMis2, 
                            NumericVector Weight,
                            NumericVector beta1, NumericVector beta2, double xi, double sigma, 
                            double gamma1, NumericVector gamma, NumericVector alpha1, NumericVector alpha0, double sigma_e){
  int i,j;
  int nsam = DesignMatrix1.nrow();
  int nval = ValidationMatrix1.nrow();
  int ncov1 = DesignMatrix1.ncol();
  int ncov2 = DesignMatrix2.ncol();
  int ncovMea = CovMis1.ncol();
  int ncovMis = CovMis2.ncol();
  NumericMatrix Vii(2,2),Viinverse(2,2);
  NumericVector eta1(nsam + nval),eta2(nsam + nval),eta2expit(nsam + nval);
  NumericVector mta_1(nsam + nval),mta_20(nsam + nval),mta_21(nsam + nval);
  NumericVector Ufunc(ncov1 + ncov2 + 2);
  NumericVector vi2(nsam + nval),Y1starstar(nsam),Y2starstar(nsam),Y1starstartilde(nsam),Y1Y2sst(nsam);
  double pi0,pi1;
  double Delta0,Delta1,Delta,S1,S2,S3;
  
 
  sigma_e = fabs(sigma_e);
  sigma = fabs(sigma);
  
  // The mismeasured data
  
  for (j = 0; j < ncov1 + ncov2 + 2; ++j){
    Ufunc(j) = 0;
  }
  
  for(i = 0; i < nsam; ++i){ 
    eta1(i) = 0;
    eta2(i) = 0;
    for (j = 0; j < ncov1; ++j){
      eta1(i) += beta1(j) * DesignMatrix1(i,j);
    }
    for (j = 0; j < ncov2; ++j){
      eta2(i) += beta2(j) * DesignMatrix2(i,j);
    }
    eta2expit(i) = exp(eta2(i))/(1+exp(eta2(i)));
    if (eta2expit(i) != eta2expit(i)) {eta2expit(i) = 1;}
    for (j = 0; j < ncovMis; ++j){
      mta_21(i) += alpha1(j) * CovMis2(i,j);
      mta_20(i) += alpha0(j) * CovMis2(i,j);
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
    
    Y2starstar(i) = (Y2star(i) - pi0) / (1 - pi0 - pi1);
    
    for (j = 0; j < ncovMea; ++j){
      if (j == 1) {
        mta_1(i) += gamma(j) * CovMis1(i,j) * Y2starstar(i);
      } else {
        mta_1(i) += gamma(j) * CovMis1(i,j);
      }
    }
    
    Y1starstar(i) = (Y1star(i) - mta_1(i)) / gamma1;
    
    
    
    Delta0 = (pi0-pow(pi0,2)) / pow(1-pi0-pi1,2);
    Delta1 = (pi1-pow(pi1,2)) / pow(1-pi0-pi1,2);
    if (Y2star(i)==1) {
      Delta = (Delta1-Delta0*pi1-Delta1*pi0) / (1-pi0-pi1);
    } else {
      Delta = (Delta0-Delta0*pi1-Delta1*pi0) / (1-pi0-pi1);
    }
    
    
    Y1starstartilde(i) =  pow(Y1starstar(i),2) - pow(sigma_e,2)/pow(gamma1,2) - pow(gamma(1),2) / pow(gamma1,2) * Delta;
    Y1Y2sst(i) = Y1starstar(i) * Y2starstar(i) + gamma(1) / gamma1 * Delta;
    
    S1 = Y1starstartilde(i) - 2 * eta1(i) * Y1starstar(i) + pow(eta1(i),2);
    S2 = Y1Y2sst(i) - eta1(i)* Y2starstar(i) - eta2expit(i)* Y1starstar(i) + eta1(i) * eta2expit(i); 
    S3 = Y2starstar(i) - 2 * eta2expit(i) * Y2starstar(i) + pow(eta2expit(i),2);
    
    for (j = 0; j < ncov1; ++j){
      Ufunc(j) += (DesignMatrix1(i,j) * Viinverse(0,0) * (Y1starstar(i) - eta1(i))
                  + DesignMatrix1(i,j) * Viinverse(0,1) * (Y2starstar(i) - eta2expit(i)))
      * Weight(j);
    }
    
    for (j = ncov1; j < ncov1 + ncov2; ++j){
      Ufunc(j) += (DesignMatrix2(i,j-ncov1) * eta2expit(i) * (1-eta2expit(i)) * Viinverse(1,0) * (Y1starstar(i) - eta1(i))
      + DesignMatrix2(i,j-ncov1) * eta2expit(i) * (1-eta2expit(i)) * Viinverse(1,1) * (Y2starstar(i) - eta2expit(i)))
      * Weight(j);
    }
    
    Ufunc(ncov1 + ncov2) += (2 * sigma * (S1-Vii(0,0)) 
      + 2*xi * pow(vi2(i),0.5) * (S2-Vii(0,1)))
      * Weight(ncov1 + ncov2);
    Ufunc(ncov1 + ncov2 + 1) += 2 * sigma * pow(vi2(i),0.5) * (S2-Vii(0,1))
      * Weight(ncov1 + ncov2 + 1);
    
  }
  
  
  // Validation Data

  for(i = nsam; i < nsam + nval; ++i){
    eta1(i) = 0;
    eta2(i) = 0;
    for (j = 0; j < ncov1; ++j){
      eta1(i) += beta1(j) * ValidationMatrix1(i-nsam,j);
    }
    for (j = 0; j < ncov2; ++j){
      eta2(i) += beta2(j) * ValidationMatrix2(i-nsam,j);
    }
    eta2expit(i) = exp(eta2(i))/(1+exp(eta2(i)));
    if (eta2expit(i) != eta2expit(i)) {eta2expit(i) = 1;}
  }

  for(i = nsam; i < nsam + nval; ++i){
    vi2(i) = eta2expit(i) * (1 - eta2expit(i));

    Vii(0,0) = pow(sigma,2);
    Vii(0,1) = xi*sigma*pow(vi2(i),0.5);
    Vii(1,0) = Vii(0,1);
    Vii(1,1) = vi2(i);

    Viinverse(0,0) = 1 / pow(sigma,2) / ( 1 - pow(xi,2));
    Viinverse(0,1) = - xi / (1 - pow(xi,2)) / sigma / pow(vi2(i),0.5);
    Viinverse(1,0) = Viinverse(0,1);
    Viinverse(1,1) = 1 / (1 - pow(xi,2)) / vi2(i);


    S1 = pow(Y1(i-nsam)-eta1(i),2);
    S2 = (Y1(i-nsam)-eta1(i))*(Y2(i-nsam)-eta2expit(i));
    S3 = pow(Y2(i-nsam)-eta2expit(i),2);

    for (j = 0; j < ncov1; ++j){
      Ufunc(j) += (ValidationMatrix1(i-nsam,j) * Viinverse(0,0) * (Y1(i-nsam) - eta1(i))
      + ValidationMatrix2(i-nsam,j) * Viinverse(0,1) * (Y2(i-nsam) - eta2expit(i))) 
      * (1- Weight(j));
    }

    for (j = ncov1; j < ncov1 + ncov2; ++j){
      Ufunc(j) += (ValidationMatrix2(i-nsam,j-ncov1) * eta2expit(i) * (1-eta2expit(i)) * Viinverse(1,0) * (Y1(i-nsam) - eta1(i))
      + ValidationMatrix2(i-nsam,j-ncov1) * eta2expit(i) * (1-eta2expit(i)) * Viinverse(1,1) * (Y2(i-nsam) - eta2expit(i)))
      * (1- Weight(j));
    }

    Ufunc(ncov1 + ncov2) += (2 * sigma * (S1-Vii(0,0)) 
      + 2*xi * pow(vi2(i),0.5) * (S2-Vii(0,1))) 
      * (1- Weight(ncov1 + ncov2));
    Ufunc(ncov1 + ncov2 + 1) += (2 * sigma * pow(vi2(i),0.5) * (S2-Vii(0,1))) 
      * (1- Weight(ncov1 + ncov2 + 1));

  }
  
  return(Ufunc);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
