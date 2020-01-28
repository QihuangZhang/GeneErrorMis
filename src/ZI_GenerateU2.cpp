# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

NumericVector ZI_GenerateU2(NumericVector Y,NumericVector U1, NumericVector mu2){
  int i;
  const int nsam = Y.length();
  
  NumericVector U2gen(nsam); 
  
  for (i = 0; i < nsam; ++ i){
    if (U1[i]==0) {U2gen[i]=R::rpois(mu2[i]);} else {U2gen[i] = Y[i];}
  }

  return U2gen; 
}

/*** R
# set.seed(2018)
# Y <- datasim$Y
# U1 <- datasim$U1
# beta2 <- c(1,0.5)
# mu2 <-  as.matrix(datasim[,c("intercept","X2")]) %*% t(t(beta2))
# 
# set.seed(2018)
# mat <- ZI_GenerateU2(Y, U1, exp(mu2))
# set.seed(2018)
# mat2 <- GenerateU2(Y, U1, exp(mu2),2000)

*/
