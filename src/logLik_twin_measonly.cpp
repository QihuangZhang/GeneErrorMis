# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double logLik_twin_measonly(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2, double rho, 
                                    NumericVector beta1, NumericVector beta2, double gamma, 
                                    double sigma, double sigma_e, double sigma_g){
  int i,j,k,t,l,s,r;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  NumericVector Y1(nsam),mta_1(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam/2);
  double eta2exp,f2i,f1i,fi,fip_p1,fip_p2;

  const double weightPi = pow(M_PI,-2.5);

  sigma = fabs(sigma);
  sigma_e = fabs(sigma_e);
  sigma_g = fabs(sigma_g);
  
  for(i = 0; i < nsam/2; ++i){
    fkkt(i)=0;
  }
  
  for(i = 0; i < nsam; ++i)
  { eta1(i)=0;
    eta2(i)=0;
    for (j = 0; j < ncov; ++j){
      eta1(i) += beta1(j) * Covariates(i,j);
      eta2(i) += beta2(j) * Covariates(i,j);
    }
  }

  for(k = 0; k < nS; ++k){
    for(l = 0; l < nS; ++l){
      for (s = 0;s < nS; ++s){
        for(i = 0; i < nsam/2; ++i)
        { // Person 1
          eta1all(i) = eta1(i) + sqrt(2*(1-rho)) * sigma_g * v(k) + sqrt(2*(rho)) * sigma_g * v(s);
          eta2exp = eta2(i) +  sqrt(2*(1-rho)) * sigma_g * v(k) + sqrt(2*(rho)) * sigma_g * v(s);
          eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
          if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}

          // Person 2
          eta1all(i+nsam/2) = eta1(i+nsam/2) + sqrt(2*(1-rho)) * sigma_g * v(l) + sqrt(2*(rho)) * sigma_g * v(s);
          eta2exp = eta2(i+nsam/2) +  sqrt(2*(1-rho)) * sigma_g * v(l) + sqrt(2*(rho)) * sigma_g * v(s);
          eta2all(i+nsam/2) = exp(eta2exp) / (1 + exp(eta2exp));
          if (eta2all(i+nsam/2) != eta2all(i+nsam/2)) { eta2all(i+nsam/2)=1 ;}
        }
        for(t = 0; t < nS; ++t){
          for (r = 0; r < nS; ++r){
            for(i = 0; i < nsam/2; ++i){
              // Person 1
              Y1(i) = eta1all(i) + sqrt(2) * sigma * v(t);

              mta_1(i)= Y1(i);
              if (isY2in == 1){
                mta_1(i) += gamma * Y2star(i);
              } else {
                mta_1(i) += gamma * (2 * Y2star(i) - 1);
              }
              f1i = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-mta_1(i))/fabs(sigma_e), 2.0 ));

              if (Y2star(i) == 1) {
                f2i = eta2all(i);
              } else {
                f2i = 1 - eta2all(i);
              }

              fip_p1 = f1i*f2i;

              // Person 2
              Y1(i+nsam/2) = eta1all(i+nsam/2) + sqrt(2) * sigma * v(r);

              mta_1(i+nsam/2)= Y1(i+nsam/2);
              if (isY2in == 1){
                mta_1(i+nsam/2) += gamma * Y2star(i+nsam/2);
              } else {
                mta_1(i+nsam/2) += gamma * (2 * Y2star(i+nsam/2) - 1);
              }
              f1i = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i+nsam/2)-mta_1(i+nsam/2))/fabs(sigma_e), 2.0 ));

              if (Y2star(i+nsam/2) == 1) {
                f2i = eta2all(i+nsam/2);
              } else {
                f2i = 1 - eta2all(i+nsam/2);
              }

              fip_p2 = f1i*f2i;

              fkkt(i) +=  w(t) * w(r) * w(k) * w(l) * w(s)   * fip_p1 * fip_p2 * weightPi;
            }
          }

        }
      }
    }

  }


  fi = 0;
  for (i = 0; i < nsam/2; ++i){
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
