# include <Rcpp.h>
using namespace Rcpp;


double gTruncPois(double lambda){
  double xstar;
  xstar = R::rpois(lambda);
  
  while (xstar==0) {
    xstar = R::rpois(lambda);
  }
  
  return xstar;
}

// [[Rcpp::export]]

NumericVector ZI_GenerateU1(NumericVector Y,NumericVector U2, NumericVector mu1){
  int i;
  const int nsam = Y.length();
  
  NumericVector U1gen(nsam); 
  
  for (i = 0; i < nsam; ++ i){
    if (Y[i]>0) {U1gen[i] = gTruncPois(mu1[i]);} else{
      if (U2[i]>0) {U1gen[i] =  0;} else {
        U1gen[i] =  R::rpois(mu1[i]);
      }
    }
  }

  return U1gen; 
}

/*** R
# set.seed(2018)
# Y <- datasim$Y
# U2 <- datasim$U2
# beta1 <- c(-0.7,0.7)
# mu1 <-  as.matrix(datasim[,c("intercept","X1")]) %*% t(t(beta1))
# 
# set.seed(2018)
# mat <- ZI_GenerateU1(Y, U2, exp(mu1))
# set.seed(2018)
# mat2 <- GenerateU1(Y, U2, exp(mu1),2000)

*/
