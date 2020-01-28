# include <Rcpp.h>
# include<iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix GEE_GAMMAInsIVI(NumericVector Y1star, NumericVector Y2star, NumericMatrix DesignMatrix1, NumericMatrix DesignMatrix2,
                             NumericMatrix CovMis1, NumericMatrix CovMis2,
                             NumericVector beta1, NumericVector beta2, double xi, double sigma,
                             double gamma1, NumericVector gamma, NumericVector alpha1, NumericVector alpha0, double sigma_e,
                             int fixgamma1, NumericVector fixgamma, int fixsigma_e, NumericVector fixalpha1, NumericVector fixalpha0){
  int i,j,k,t;
  int nsam = DesignMatrix1.nrow();
  int ncov1 = DesignMatrix1.ncol();
  int ncov2 = DesignMatrix2.ncol();
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
  
  NumericMatrix Vii(2,2),Viinverse(2,2);
  NumericMatrix Dmatrix1(ncov1+ncov2+2,ncov1+ncov2+2),Dmatrix2(ncov1+ncov2+2,ncov1+ncov2+2+nq);
  NumericMatrix GammaMatrixI(ncov1+ncov2+2+nq,ncov1+ncov2+2+nq);
  NumericVector eta1(nsam),eta2(nsam),eta2expit(nsam);
  NumericVector mta_1(nsam),mta_20(nsam),mta_21(nsam);
  NumericVector vi2(nsam),Y1starstar(nsam),Y2starstar(nsam),Y1starstartilde(nsam),Y1Y2sst(nsam);
  double pi0,pi1;
  double Delta0,Delta1,Delta;
  double Del1,Del2,Del3,Del4;
  int index;
  
  sigma = fabs(sigma);
  
  // Evaluate GammaMatrixI
  for (j = 0; j < ncov1 + ncov2 + 2 + nq; ++j){
    for (k = 0; k < ncov1 + ncov2 + 2 + nq; ++k){
      GammaMatrixI(j,k) = 0;
    }
  }
  
  for (j = 0; j < ncov1 + ncov2 + 2; ++j){
    for (k = 0; k < ncov1 + ncov2 + 2; ++k){
      Dmatrix1(j,k) = 0;
      Dmatrix2(j,k) = 0;
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
    
    for(j = 0; j < ncov1; ++j){
      Dmatrix1(j,0) = DesignMatrix1(i,j) * Viinverse(0,0);
      Dmatrix1(j,1) = DesignMatrix1(i,j) * Viinverse(0,1);
    }
    for(j = ncov1; j < ncov1+ncov2; ++j){
      Dmatrix1(j,0) = DesignMatrix2(i,j-ncov1) * vi2(i) * Viinverse(1,0);
      Dmatrix1(j,1) = DesignMatrix2(i,j-ncov1) * vi2(i) * Viinverse(1,1);
    }
    
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
    
    // Construct the Dmatrix2 additional part of mismeasurement parameters
    index = ncov1+ncov2+1;
    
    if (fixgamma1!=1){
      index = index + 1;
      Del1 = - Y1starstar(i)/gamma1;
      Del2 = 0;
      Del3 = -2/gamma1 * (pow(Y1starstartilde(i),2)-pow(sigma_e/gamma1,2)-pow(gamma(1)/gamma1,2)*Delta);
      Del4 = - Y1Y2sst(i)/gamma1;
      Dmatrix2(0,index) = Del1;
      Dmatrix2(1,index) = Del2;
      Dmatrix2(2,index) = Del3 - 2 * eta1(i) * Del1;
      Dmatrix2(3,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
      Dmatrix2(4,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
      Dmatrix2(5,index) = (1-2*eta2expit(i)) * Del2;
    }
    
    for(j = 0; j < gamma.length(); ++j){
      if (fixgamma(j)!=1){
        index = index + 1;
        if (j == 1) {
          Del1 = - Y2starstar(i)/gamma1;
          Del2 = 0;
          Del3 = -2/gamma1 * Y1Y2sst(i);
          Del4 = -pow(Y2starstar(i),2)/gamma1 + Delta / gamma1;
        } else {
          Del1 = - CovMis1(i,j)/gamma1;
          Del2 = 0;
          Del3 = -2/gamma1 * Y1starstar(i) * CovMis1(i,j);
          Del4 = - Y2starstar(i) / gamma1 * CovMis1(i,j);
        }
        
        Dmatrix2(0,index) = Del1;
        Dmatrix2(1,index) = Del2;
        Dmatrix2(2,index) = Del3 - 2 * eta1(i) * Del1;
        Dmatrix2(3,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
        Dmatrix2(4,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
        Dmatrix2(5,index) = (1-2*eta2expit(i)) * Del2;
      }
    } 
    
    if (fixsigma_e!=1){
      index = index + 1;
      Del1 = 0;
      Del2 = 0;
      Del3 = -2*sigma_e/pow(gamma1,2);
      Del4 = 0;
      Dmatrix2(0,index) = Del1;
      Dmatrix2(1,index) = Del2;
      Dmatrix2(2,index) = Del3 - 2 * eta1(i) * Del1;
      Dmatrix2(3,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
      Dmatrix2(4,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
      Dmatrix2(5,index) = (1-2*eta2expit(i)) * Del2;
    }
    
    for(j = 0; j < fixalpha1.length(); ++j){
      if (fixalpha1(j)!=1) {
        index = index + 1;
        
        Del1 = gamma(1)/gamma1*(-Y2starstar(i))/(1-pi0-pi1)*(CovMis2(i,j)*exp(mta_21(i)))/pow(1+exp(mta_21(i)),2);
        Del2 = Y2starstar(i)/(1-pi0-pi1)*(CovMis2(i,j)*exp(mta_21(i)))/pow(1+exp(mta_21(i)),2);
        Del3 = (-0.2e1 * Y1starstar(i) * (Y2star(i) - pi0) * gamma(1) / gamma1 * pow((1 - pi0 - pi1), (-2)) - gamma(1) * gamma(1) * pow(gamma1, -0.2e1) * (pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * (1 - Y2star(i)) * ((-2 * pi1 + 1) * pow((1 - pi0 - pi1), (-2)) + 2 * (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-3))) / (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), 2) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) + 2 * pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) * Y2star(i) / (1 - pi0 - pi1) - 2 * (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-3)) * pi1 - (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2)) - (-2 * pi1 + 1) * pow((1 - pi0 - pi1), (-2)) * pi0 - 2 * (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-3)) * pi0) / (1 - pi0 - pi1) - gamma(1) * gamma(1) * pow(gamma1, -0.2e1) * (pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) - (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2)) * pi1 - (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2)) * pi0) * pow((1 - pi0 - pi1), (-2))
        )*(CovMis2(i,j)*exp(mta_21(i)))/pow(1+exp(mta_21(i)),2);
        Del4 = (-pow(Y2star(i) - pi0, 0.2e1) * pow(0.1e1 - pi0 - pi1, -0.3e1) * gamma(1) / gamma1 + Y1starstar(i) * (Y2star(i) - pi0) * pow(0.1e1 - pi0 - pi1, -0.2e1) + gamma(1) / gamma1 * (pow((-pi1 * pi1 + pi1) * pow(0.1e1 - pi0 - pi1, -0.2e1), 0.1e1 - Y2star(i)) * (0.1e1 - Y2star(i)) * ((-0.2e1 * pi1 + 0.1e1) * pow(0.1e1 - pi0 - pi1, -0.2e1) + 0.2e1 * (-pi1 * pi1 + pi1) * pow(0.1e1 - pi0 - pi1, -0.3e1)) / (-pi1 * pi1 + pi1) * pow(0.1e1 - pi0 - pi1, 0.2e1) * pow((-pi0 * pi0 + pi0) * pow(0.1e1 - pi0 - pi1, -0.2e1), Y2star(i)) + 0.2e1 * pow((-pi1 * pi1 + pi1) * pow(0.1e1 - pi0 - pi1, -0.2e1), 0.1e1 - Y2star(i)) * pow((-pi0 * pi0 + pi0) * pow(0.1e1 - pi0 - pi1, -0.2e1), Y2star(i)) * Y2star(i) / (0.1e1 - pi0 - pi1) - 0.2e1 * (-pi0 * pi0 + pi0) * pow(0.1e1 - pi0 - pi1, -0.3e1) * pi1 - (-pi0 * pi0 + pi0) * pow(0.1e1 - pi0 - pi1, -0.2e1) - (-0.2e1 * pi1 + 0.1e1) * pow(0.1e1 - pi0 - pi1, -0.2e1) * pi0 - 0.2e1 * (-pi1 * pi1 + pi1) * pow(0.1e1 - pi0 - pi1, -0.3e1) * pi0) / (0.1e1 - pi0 - pi1) + gamma(1) / gamma1 * (pow((-pi1 * pi1 + pi1) * pow(0.1e1 - pi0 - pi1, -0.2e1), 0.1e1 - Y2star(i)) * pow((-pi0 * pi0 + pi0) * pow(0.1e1 - pi0 - pi1, -0.2e1), Y2star(i)) - (-pi0 * pi0 + pi0) * pow(0.1e1 - pi0 - pi1, -0.2e1) * pi1 - (-pi1 * pi1 + pi1) * pow(0.1e1 - pi0 - pi1, -0.2e1) * pi0) * pow(0.1e1 - pi0 - pi1, -0.2e1)
        )*(CovMis2(i,j)*exp(mta_21(i)))/pow(1+exp(mta_21(i)),2);
        
        Dmatrix2(0,index) = Del1;
        Dmatrix2(1,index) = Del2;
        Dmatrix2(2,index) = Del3 - 2 * eta1(i) * Del1;
        Dmatrix2(3,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
        Dmatrix2(4,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
        Dmatrix2(5,index) = (1-2*eta2expit(i)) * Del2;
      }
    } 
    
    for(j = 0; j < fixalpha0.length(); ++j){
        if (fixalpha0(j)!=1) {
          index = index +1;
          
          Del1 = gamma(1)/gamma1*(1-Y2starstar(i))/(1-pi0-pi1)*(CovMis2(i,j)*exp(mta_20(i)))/pow(1+exp(mta_20(i)),2);
          Del2 = (Y2starstar(i)-1)/(1-pi0-pi1)*(CovMis2(i,j)*exp(mta_20(i)))/pow(1+exp(mta_20(i)),2);
          Del3 =(0.2e1 * Y1starstar(i) * (0.1e1 / (1 - pi0 - pi1) * gamma(1) - (Y2star(i) - pi0) * pow((1 - pi0 - pi1), (-2)) * gamma(1)) / gamma1 - gamma(1) * gamma(1) * pow(gamma1, -0.2e1) * (2 * pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * (1 - Y2star(i)) / (1 - pi0 - pi1) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) + pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) * Y2star(i) * ((-2 * pi0 + 1) * pow((1 - pi0 - pi1), (-2)) + 2 * (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-3))) / (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), 2) - (-2 * pi0 + 1) * pow((1 - pi0 - pi1), (-2)) * pi1 - 2 * (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-3)) * pi1 - 2 * (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-3)) * pi0 - (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))) / (1 - pi0 - pi1) - gamma(1) * gamma(1) * pow(gamma1, -0.2e1) * (pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) - (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2)) * pi1 - (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2)) * pi0) * pow((1 - pi0 - pi1), (-2))
          )*(CovMis2(i,j)*exp(mta_20(i)))/pow(1+exp(mta_20(i)),2);
          Del4 = ((0.1e1 / (1 - pi0 - pi1) * gamma(1) - (Y2star(i) - pi0) * pow((1 - pi0 - pi1), (-2)) * gamma(1)) / gamma1 * (Y2star(i) - pi0) / (1 - pi0 - pi1) - (Y1starstar(i) / (1 - pi0 - pi1)) + (Y1starstar(i) * (Y2star(i) - pi0) * pow((1 - pi0 - pi1), (-2))) + gamma(1) / gamma1 * (2 * pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * (1 - Y2star(i)) / (1 - pi0 - pi1) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) + pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) * Y2star(i) * ((-2 * pi0 + 1) * pow((1 - pi0 - pi1), (-2)) + 2 * (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-3))) / (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), 2) - (-2 * pi0 + 1) * pow((1 - pi0 - pi1), (-2)) * pi1 - 2 * (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-3)) * pi1 - 2 * (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-3)) * pi0 - (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))) / (1 - pi0 - pi1) + gamma(1) / gamma1 * (pow(((-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2))), (1 - Y2star(i))) * pow(((-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2))), Y2star(i)) - (-pi0 * pi0 + pi0) * pow((1 - pi0 - pi1), (-2)) * pi1 - (-pi1 * pi1 + pi1) * pow((1 - pi0 - pi1), (-2)) * pi0) * pow((1 - pi0 - pi1), (-2))
          )*(CovMis2(i,j)*exp(mta_20(i)))/pow(1+exp(mta_20(i)),2);
            
            Dmatrix2(0,index) = Del1;
            Dmatrix2(1,index) = Del2;
            Dmatrix2(2,index) = Del3 - 2 * eta1(i) * Del1;
            Dmatrix2(3,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
            Dmatrix2(4,index) = Del4 - eta1(i) * Del2 - eta2expit(i) * Del1;
            Dmatrix2(5,index) = (1-2*eta2expit(i)) * Del2;
        }
      } 
    
    
    for(j = 0; j < ncov1 + ncov2 + 2; ++j)
    {
      for(k = 0; k < ncov1 + ncov2 + 2 + nq; ++k)
      {
        for(t=0; t<  ncov1 + ncov2 + 2; ++t)
        {
          GammaMatrixI(j,k) += Dmatrix1(j,t) * Dmatrix2(t,k);
        }
      }
    }
    
  }
  
  return(GammaMatrixI);
}
                                         
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
