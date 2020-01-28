# include <Rcpp.h>
# include<iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix GEE_GAMMAIns(NumericVector Y1star, NumericVector Y2star, NumericMatrix DesignMatrix1, NumericMatrix DesignMatrix2,
                                    NumericVector beta1, NumericVector beta2, double xi, 
                                    double sigma){
  int i,j,k,t;
  int nsam = DesignMatrix1.nrow();
  int ncov1 = DesignMatrix1.ncol();
  int ncov2 = DesignMatrix2.ncol();
  NumericMatrix Vii(2,2),Viinverse(2,2);
  NumericMatrix Dmatrix1(ncov1+ncov2+2,ncov1+ncov2+2),Dmatrix2(ncov1+ncov2+2,ncov1+ncov2+2);
  NumericMatrix GammaMatrix(ncov1+ncov2+2,ncov1+ncov2+2);
  NumericVector eta1(nsam),eta2(nsam),eta2expit(nsam);
  NumericVector mta_1(nsam),mta_20(nsam),mta_21(nsam);
  NumericVector vi2(nsam),Y1starstar(nsam),Y2starstar(nsam),Y1starstartilde(nsam),Y1Y2sst(nsam);
  
  // Delete this after debugging
  NumericVector temp(nsam);
  
  sigma = fabs(sigma);
  
  for (j = 0; j < ncov1 + ncov2 + 2; ++j){
    for (k = 0; k < ncov1 + ncov2 + 2; ++k){
      Dmatrix1(j,k) = 0;
      Dmatrix2(j,k) = 0;
      GammaMatrix(j,k) = 0;
    }
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
    
    for(j = 0; j < ncov1; ++j){
      Dmatrix1(j,0) = DesignMatrix1(i,j) * Viinverse(0,0);
      Dmatrix1(j,1) = DesignMatrix1(i,j) * Viinverse(0,1);
    }
    for(j = ncov1; j < ncov1+ncov2; ++j){
      Dmatrix1(j,0) = DesignMatrix2(i,j-ncov1) * vi2(i) * Viinverse(1,0);
      Dmatrix1(j,1) = DesignMatrix2(i,j-ncov1) * vi2(i) * Viinverse(1,1);
    }
    
    Dmatrix1(ncov1+ncov2,2) = 2 * sigma;
    Dmatrix1(ncov1+ncov2,3) = xi * pow(vi2(i),0.5);
    Dmatrix1(ncov1+ncov2,4) = xi * pow(vi2(i),0.5);
    
    Dmatrix1(ncov1+ncov2+1,2) = 0;
    Dmatrix1(ncov1+ncov2+1,3) = sigma  * pow(vi2(i),0.5);
    Dmatrix1(ncov1+ncov2+1,4) = sigma  * pow(vi2(i),0.5);
    Dmatrix1(ncov1+ncov2+1,5) = 0;
    
    // Construct the Dmatrix2
    for(j = 0; j < ncov1; ++j){
      Dmatrix2(0,j) = - DesignMatrix1(i,j);
    }
    
    for(j = ncov1; j < ncov1+ncov2; ++j){
      Dmatrix2(1,j) = - DesignMatrix2(i,j-ncov1) * vi2(i);
    }
    
    for(j = ncov1; j < ncov1 + ncov2; ++j){
      Dmatrix2(3,j) = - 0.5 * xi * sigma * (1-2 * eta2expit(i)) * DesignMatrix2(i,j-ncov1) * pow(vi2(i),0.5);
      Dmatrix2(4,j) = Dmatrix2(3,j);
      Dmatrix2(5,j) = - (1-2 * eta2expit(i)) * DesignMatrix2(i,j-ncov1) * vi2(i);
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
          GammaMatrix(j,k) += Dmatrix1(j,t) * Dmatrix2(t,k) ;
        }
      }
    }
    
  }
  
  return(GammaMatrix);
}
                                         
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
