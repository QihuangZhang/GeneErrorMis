# include <Rcpp.h>
using namespace Rcpp;

double ajumpfunc(double theta, double sigma){
  return R::rnorm(theta,sigma);
}


// This function based on the assumption that Zminus follows a logsitic regression


// [[Rcpp::export]]

NumericVector ZI_GenerateAlphaMNMetro(NumericVector Par, NumericVector Y, NumericVector Z, NumericMatrix Covar, NumericVector propsigma, 
                                       NumericVector priormu, NumericVector priorSigmas){
  int i,j,k;
  int nsam=Covar.nrow();
  int ncov=Covar.ncol();
  
  double u, parastar, logratioposterior;
  
  NumericVector Pargen(ncov);
  NumericVector mu(nsam), muprop(nsam), pi(nsam), piprop(nsam);
  
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
    parastar = ajumpfunc(Pargen[i], propsigma[i]);
    
    logratioposterior = (Pargen[i]-priormu[i])*(Pargen[i]-priormu[i])/(2*priorSigmas[i]*priorSigmas[i])
      -(parastar-priormu[i])*(parastar-priormu[i])/(2*priorSigmas[i]*priorSigmas[i]);
    
    // std::cout << 'P' << 'r' << logratioposterior  << ' ';  
      
    for (j = 0; j < nsam; ++j){
      muprop[j] = mu[j] + Covar(j,i) * (parastar - Pargen[i]);
      // pi[j] = exp(mu[j])/(1+exp(mu[j]));
      // piprop[j] = exp(muprop[j])/(1+exp(muprop[j]));
      // logratioposterior += Z[j]*log(piprop[j]) + (Y[j]-Z[j])*log(1-piprop[j])
      //   - Z[j]*log(pi[j]) - (Y[j]-Z[j])*log(1-pi[j]);
      logratioposterior += Z[j] * (muprop[j] - mu[j]) - Y[j] * (log(1 + exp(muprop[j]))
                                  - log(1 + exp(mu[j])));
    }
    
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

# alphaminus0 <- -0.5
# Y <- datasim$Y
# Zminus <- datasim$Zminus
# Covarminus <- datasim$Xminus
# priormu <-0
# 
# priorSigma <- 10000

# alphaminus0 <-0.1
# # record1 <- NULL
# # for (i in 1:50000){
# #   alphaminus <- GenerateAlphaminusCpp(alphaminus, Y, as.matrix(Covarminus), Zminus, priormu, priorSigma)
# #   record1 <- c(record1,alphaminus)
# # }
# 
# alphaminus <- alphaminus0
# record2<- NULL
# for (i in 1:50000){
#   alphaminus <-  ZI_GenerateAlphaMNMetro(alphaminus, Y, Zminus, as.matrix(Covarminus),  0.05,
#                                          priormu,  priorSigmas=priorSigma)
#   record2<- c(record2, alphaminus)
# }
#
# # hist(record1)
# hist(record2[1000:50000])

# alphaminus <- c(0,0)
# Y <- 
# 
# 
# alphaminus <- ZI_GenerateAlphaMNMetro(alphaminus, Y, Zminus, as.matrix(Covarminus),  propsigma= propsigmai[7:8],  priormu,  priorSigmas=priorSigma)
*/
