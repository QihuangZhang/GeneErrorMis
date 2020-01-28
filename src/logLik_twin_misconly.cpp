# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double logLik_twin_misconly(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2, double rho, 
                                    NumericVector beta1, NumericVector beta2, NumericVector alpha1, NumericVector alpha0,
                                    double sigma, double sigma_g){
  int i,j,k,l,s;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  int ncovMis=CovMis2.ncol();
  NumericVector mta_20(nsam),mta_21(nsam),mta_20expit(nsam),mta_21expit(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam/2);
  double eta2exp,fip_1p,fip_2p,f2i1,f2i0,f1i,fi;
  
  const double weightPi = pow(M_PI,-0.5*3);
  
  sigma = fabs(sigma);
  sigma_g = fabs(sigma_g);
  
  for(i = 0; i < nsam/2; ++i){
    fkkt(i)=0;
  }
  
  for(i = 0; i < nsam; ++i)
  { eta1(i)=0;
    eta2(i)=0;
    mta_21(i) = 0;
    mta_20(i) = 0;
    for (j = 0; j < ncov; ++j){
      eta1(i) += beta1(j) * Covariates(i,j);
      eta2(i) += beta2(j) * Covariates(i,j);
    }
    for (j = 0; j < ncovMis; ++j){
      mta_21(i) += alpha1(j) * CovMis2(i,j);
      mta_20(i) += alpha0(j) * CovMis2(i,j);
    }
    mta_20expit(i) = exp(mta_20(i))/(1+exp(mta_20(i)));
    mta_21expit(i) = exp(mta_21(i))/(1+exp(mta_21(i)));
    if (mta_20expit(i) != mta_20expit(i)) {mta_20expit(i) = 1;}
    if (mta_21expit(i) != mta_21expit(i)) {mta_21expit(i) = 1;}
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
        for(i = 0; i < nsam/2; ++i){
          // Person 1
          f1i = ( 1 / (fabs(sigma) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-eta1all(i))/fabs(sigma), 2.0 ));
          
          f2i1 = pow(mta_21expit(i), 1-Y2star(i)) * pow(1-mta_21expit(i), Y2star(i)) ;
          f2i0 = pow(mta_20expit(i), Y2star(i)) * pow(1-mta_20expit(i), 1-Y2star(i)) ;
          
          fip_1p = f1i * (f2i0 * (1-eta2all(i)) + f2i1 * eta2all(i));
          
          // Person 2
          f1i = ( 1 / (fabs(sigma) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i+nsam/2)-eta1all(i+nsam/2))/fabs(sigma), 2.0 ));
          
          f2i1 = pow(mta_21expit(i+nsam/2), 1-Y2star(i+nsam/2)) * pow(1-mta_21expit(i+nsam/2), Y2star(i+nsam/2)) ;
          f2i0 = pow(mta_20expit(i+nsam/2), Y2star(i+nsam/2)) * pow(1-mta_20expit(i+nsam/2), 1-Y2star(i+nsam/2)) ;
          
          fip_2p = f1i * (f2i0 * (1-eta2all(i+nsam/2)) + f2i1 * eta2all(i+nsam/2));
          
          fkkt(i) +=  w(k) * w(l) * w(s) * fip_1p * fip_2p * weightPi;
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
