# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

NumericVector ZI_GeneratePoiPar_Binary(NumericVector Par, NumericMatrix covariates, NumericVector Outcome, 
                                       NumericVector priorgamma){
  int i,j,k;
  int nsam=covariates.nrow();
  int ncov=covariates.ncol();
  
  double rate, shape;
  
  NumericVector Pargen(ncov), covsumtilde(nsam);
  
  

  for (i = 0; i < ncov; ++ i){
    for (k = 0; k < nsam; ++ k){
      covsumtilde[k] = 0;
    }
    
    rate = 0; shape = 0;
    for (j = 0; j < ncov; ++ j){
      if (j != i) {
        for (k = 0; k < nsam; ++ k){
          covsumtilde[k] += covariates(k,j) * Par[j];
        }
      }
    }
    
    for (k = 0; k < nsam; ++ k){
      if (covariates(k,i) == 1) {
        shape += Outcome[k];
        rate +=  exp(covsumtilde[k]);
      }
    }
    
    shape +=  priorgamma[0];
    rate += priorgamma[1];
    
    Pargen[i]  = log ( R::rgamma(shape,1/rate));
  }

  return Pargen;
}


/*** R
# set.seed(2018)
# ZI_GeneratePoiPar_Binary(c(-1,0.5), as.matrix(datasim[,c("intercept","Xplus")]),datasim$Y,c(0.01,0.01))
*/
