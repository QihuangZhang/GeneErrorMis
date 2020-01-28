# include <Rcpp.h>
using namespace Rcpp;


double gTruncRobert(double trunca){
  double u,z,rho,alphastar;
  
  alphastar = (trunca + sqrt(pow(trunca,2) + 4))/2;
  z = R::rexp(1/alphastar) + trunca;
  rho = exp(-pow(z-alphastar,2)/2);
  u = R::runif(0,1);
  while(u > rho){
    z =  R::rexp(1/alphastar) + trunca;
    rho = exp(-pow(z-alphastar,2)/2);
    u = R::runif(0,1);
  }
  return z;
}

double gTruncNorm(double a,double mu,int sign){
  double u;
  
  if (sign>0) {
    if (mu >= a) {
      u = R::rnorm(mu,1);
      while (u < a ) {
        u = R::rnorm(mu,1);
      }
      return u;
    } else{
      double Pos = gTruncRobert(a-mu);
      return Pos+mu;
    }
  } else{
    double Neg = gTruncNorm(-a,-mu,-sign);
    return -Neg;
  }
}



// [[Rcpp::export]]

NumericVector ZI_GenerateV(NumericVector Y,NumericVector Zminus, NumericVector muZminus, NumericMatrix covminus){
  
  int i,j,k;
  const int nsam = Y.length();
  const int ncov = covminus.ncol();
  int Ysum=0, index=0;
  
  for (i = 0; i < nsam; ++ i){
    Ysum += Y[i];
  }
  
  NumericMatrix Vmat(Ysum,ncov+1);
  
  for (i = 0; i < nsam; ++ i){
    if (Y[i]>0){
      if (Zminus[i]==0) {
        for (j = 0; j < Y[i]; ++ j){
          Vmat(index,0) = gTruncNorm(0,muZminus[i],-1);
          for (k = 1; k < ncov+1; ++ k){
            Vmat(index,k) = covminus(i,k-1);
          }
          ++index;
        }
      } else {
        if (Zminus[i] == Y[i]){
          for (j = 0; j < Y[i]; ++ j){
            Vmat(index,0) =  gTruncNorm(0,muZminus[i],1);
            for (k = 1; k < ncov+1; ++ k){
              Vmat(index,k) =  covminus(i,k-1);
            }
            ++index;
          }
        }else {
          for (j = 0;  j < Zminus[i]; ++j){
            Vmat(index,0) =  gTruncNorm(0,muZminus[i],1);
            for (k = 1; k < ncov+1; ++ k){
              Vmat(index,k) = covminus(i,k-1);
            }
            ++index;
          }
          for (j = Zminus[i]; j < Y[i]; ++j){
            Vmat(index,0) = gTruncNorm(0,muZminus[i],-1);
            for (k = 1; k < ncov+1; ++ k){
              Vmat(index,k) = covminus(i,k-1);
            }
            ++index;
          }
        }
      }
    }
  }
  

  return Vmat; 
}

/*** R
# set.seed(2018)
# Y <- datasim$Y
# Zminus <- datasim$Zminus
# muZminus <- datasim$intercept * (-1)
# covminus <- as.matrix(datasim$intercept)
# 
# mat <- GenerateV(Y, Zminus, muZminus, covminus)
*/
