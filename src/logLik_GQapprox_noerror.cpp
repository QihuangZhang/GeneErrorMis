# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double logLik_GQapprox_noerror(NumericVector w, NumericVector v, NumericVector Y1, NumericVector Y2, NumericMatrix Covariates, NumericVector R, 
                                    NumericVector beta1, NumericVector beta2,
                                    double sigma, double sigma_g){
  int i,j,k;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam);
  double eta2exp,fip,f1i,f2i,fi;
  
  sigma = fabs(sigma);
  sigma_g = fabs(sigma_g);
  
  for(i = 0; i < nsam; ++i)
  { fkkt(i)=0;
    eta1(i)=0;
    eta2(i)=0;
    for (j = 0; j < ncov; ++j){
      eta1(i) += beta1(j) * Covariates(i,j);
      eta2(i) += beta2(j) * Covariates(i,j);
    }
  }

  for(k = 0; k < nS; ++k){
    for(i = 0; i < nsam; ++i)
    { eta1all(i) = eta1(i) + sqrt(2*R(i)) * sigma_g * v(k);
      eta2exp = eta2(i) +  sqrt(2*R(i)) * sigma_g * v(k);
      eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
      if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}
    }
    for(i = 0; i < nsam; ++i){
        f1i = ( 1 / (fabs(sigma) * sqrt(2*M_PI))) * exp( -0.5 * pow( (Y1(i)-eta1all(i))/fabs(sigma), 2.0 ));
        f2i = pow(eta2all(i), Y2(i)) * pow(1-eta2all(i), 1-Y2(i)) ;
              
        fip = f1i * f2i;

        fkkt(i) +=  w(k)*fip / M_PI;
      }
    
  }


  fi = 0;
  for (i = 0; i < nsam; ++i){
    fi += log(fkkt(i));
  }
  
  return(fi);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# 
# sigma_g<-2
# set.seed(2018)
# u<-mvrnorm(n=1000,mu=rep(0,dim(RR_new)[1]),Sigma=RR_new*(sigma_g^2))
# random1<-matrix(rnorm(nsample*nrepeat,mean=0,sd=abs(sigma)),ncol=nrepeat)
# random2<-matrix(runif(nsample*nrepeat,min=0,max=1),ncol=nrepeat)
# 
# nsam<-dim(Covariates)[1]
# 
# beta1<-c(1,1,1)
# beta2<-c(1,1,0.1)
#   lambda<-c(1,1)
#   gamma<-c(14,0.1,0.1)
#   alpha0<-c(0.01,0.01)
#   alpha1<-c(0.01,0.01)
#   sigma<-sigma_e<-sigma_g<-2
# 
# CovMea1<- as.matrix(cbind(Covariates_all[,c("Intercept")]))
# CovMis2<- as.matrix(cbind(Covariates_all[,c("Intercept","BMD")]))
# Covariates<- as.matrix(cbind(Covariates_all[,c("Intercept","X","BW")])
# 
# Re<-logLikelihood_indep_c(u,random1,random2,Y1star,Y2star,Covariates,CovMea1,1,1,CovMis2,
#                      beta1,beta2,lambda,gamma,alpha1,alpha0,sigma,sigma_e)
# Re
# set.seed(2018)
*/
