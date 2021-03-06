# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double logLikelihood_indep_nullmeas_c(NumericMatrix u, NumericMatrix random1, NumericMatrix random2, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2,  
                                    NumericVector beta1, NumericVector beta2, NumericVector lambda, double gamma, NumericVector alpha1, NumericVector alpha0, double sigma, double sigma_e){
  int i,j,k,t;
  int nk=u.nrow();
  int nrepeat=random1.ncol();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  int ncovMis=CovMis2.ncol();
  NumericVector Y1(nsam),Y2(nsam),mta_1(nsam),debug(nsam);
  NumericVector mta_20(nsam),mta_21(nsam),mta_20expit(nsam),mta_21expit(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam),count(nsam);
  double eta2exp,f1i,f2i,fip,fi;
  
  sigma = fabs(sigma);
  
  for(i = 0; i < nsam; ++i)
  { fkkt(i)=0;
    count(i)=0;
    eta1(i)=0;
    eta2(i)=0;
    mta_21(i) = 0;
    mta_20(i) = 0;
    for (j = 0; j < ncov; ++j){
      eta1(i) += beta1(j)*Covariates(i,j);
      eta2(i) += beta2(j)*Covariates(i,j);
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
  
  for(k = 0; k < nk; ++k){
    for(i = 0; i < nsam; ++i)
    { eta1all(i) = eta1(i) + lambda(0) * u(k,i);
      eta2exp = eta2(i) + lambda(1) * u(k,i);
      eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
      if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}
    }
    for(t = 0; t < nrepeat; ++t){
      for(i = 0; i < nsam; ++i){
        Y1(i) = eta1all(i) + random1(i,t);
        if (eta2all(i) >= random2(i,t)) {Y2(i) = 1;}
        else {Y2(i) = 0;}
        
        mta_1(i)= Y1(i);
        if (isY2in == 1){
          mta_1(i) += gamma * Y2(i);
        } else {
          mta_1(i) += gamma * (2 * Y2(i) - 1);
        }
        
        if (Y2(i)==1) {f2i = pow(mta_21expit(i), 1-Y2star(i)) * pow(1-mta_21expit(i), Y2star(i)) ;}
        else {f2i = pow(mta_20expit(i), Y2star(i)) * pow(1-mta_20expit(i), 1-Y2star(i)) ;}
        f1i = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-mta_1(i))/fabs(sigma_e), 2.0 ));
        fip = f1i * f2i;
        
        if (fip == fip) {
          count(i) +=1;
          fkkt(i) +=  fip;}
      }
    }
  }
  
  
  fi = 0;
  for (i = 0; i < nsam; ++i){
    fi += log(fkkt(i)/count(i));
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
