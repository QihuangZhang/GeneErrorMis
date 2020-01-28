# include <Rcpp.h>
using namespace Rcpp;



double jumpfunc(double theta, double sigma){
  return R::rnorm(theta,sigma);
}



// [[Rcpp::export]]

NumericVector ZI_GenerateBetaMetro(NumericVector Par, NumericVector U, NumericMatrix Covar, NumericVector propsigma, 
                                       NumericVector priorgamma){
  int i,j,k;
  int nsam=Covar.nrow();
  int ncov=Covar.ncol();
  
  double u, parastar, Ucovar, expsum, logratioposterior;
  
  NumericVector Pargen(ncov);
  NumericVector mu(nsam), muprop(nsam);
  
  for (i = 0; i < ncov; ++ i){
    Pargen[i] = Par[i];
  }
  
  for (i = 0; i < nsam; ++ i){
    mu[i] = 0;
    for (k = 0; k < ncov; ++ k){
      mu[i] += Covar(i,k) * Pargen[k];
    } 
  }
  
  for (i = 0; i < ncov; ++ i){
    parastar = jumpfunc(Pargen[i], propsigma[i]);
    
    for (j = 0; j < nsam; ++j){
      muprop[j] = mu[j] + Covar(j,i) * (parastar - Pargen[i]);
    }
    
    Ucovar = 0;
    expsum = 0;
    for (j = 0; j < nsam; ++j){
      Ucovar += U[j] * Covar(j,i);
      expsum += exp(mu[j])-exp(muprop[j]);
    }
    
    logratioposterior = (Ucovar + priorgamma[0]-1) * (parastar - Pargen[i]) +
      expsum + priorgamma[1] * (exp(Pargen[i])-exp(parastar));
    
    // std::cout << 'P' << logratioposterior  << ' ';
    
    if (logratioposterior>=0) {
      Pargen[i] = parastar;
      for (k = 0; k < nsam; ++ k){
        mu[k] = muprop[k];
      }
    } else{
      u = R::runif(0,1); 
      if (log(u)<=logratioposterior) {
        Pargen[i] = parastar;
        // mu = muprop;
        for (k = 0; k < nsam; ++ k){
          mu[k] = muprop[k];
        }
      }
    }
    
  }

  return Pargen;
}


/*** R
# set.seed(2018)
# Pare <- c(0,0.5)
# U <- datasim$U2
# Covar <- datasim[,c("intercept","X2")]
# propsigma <- c(0.0001,0.0001)
# priorgamma <- c(0.001,0.001)
# 
# set.seed(2)
#  GenerateBetaMetro(Pare, U, Covar, propsigma, priorgamma)
#  Pare
#  GenerateBetaMetro(Pare, U, Covar, propsigma, priorgamma)
# 
# set.seed(2)
# ZI_GenerateBetaMetro(Pare, U, as.matrix(Covar), propsigma, priorgamma)
# Pare
# ZI_GenerateBetaMetro(Pare, U, as.matrix(Covar), propsigma, priorgamma)
*/
