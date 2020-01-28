# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector IL_score_misconly(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2, NumericVector R,  
                                    NumericVector beta1, NumericVector beta2, NumericVector alpha1, NumericVector alpha0,
                                    double sigma, double sigma_g){
  int i,j,k,t;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  int ncovMis=CovMis2.ncol();
  NumericVector mta_20(nsam),mta_21(nsam),mta_20expit(nsam),mta_21expit(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam);
  double eta2exp,uk;
  double fip,fip0,fip1,f2i1,f2i0,f1i,weightjk;
  
  const int indexbeta1 = beta1.length();
  const int indexbeta2 = indexbeta1+beta2.length();
  const int indexalpha1 = indexbeta2+alpha1.length();
  const int indexalpha0 = indexalpha1+alpha0.length();
  const int indexsigma = indexalpha0+1;
  const int indexsigmag = indexsigma+1;
  
  NumericMatrix scores(nsam,indexsigmag);
  NumericVector scoresum(indexsigmag);
  
  sigma = fabs(sigma);
  sigma_g = fabs(sigma_g);
  
  for(i = 0; i < nsam; ++i)
  { fkkt(i)=0;
    eta1(i)=0;
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
  
  for(i = 0; i < nsam; ++i){
    for(k = 0; k < nS; ++k){
      uk = sqrt(2*R(i)) * sigma_g * v(k);
      eta1all(i) = eta1(i) +  uk;
      eta2exp = eta2(i) +  uk;
      eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
      if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}
      
      f1i = ( 1 / (fabs(sigma) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-eta1all(i))/fabs(sigma), 2.0 ));
        
      f2i1 = pow(mta_21expit(i), 1-Y2star(i)) * pow(1-mta_21expit(i), Y2star(i)) ;
      f2i0 = pow(mta_20expit(i), Y2star(i)) * pow(1-mta_20expit(i), 1-Y2star(i)) ;
        
      fip0 = f1i * f2i0 * (1-eta2all(i));
      fip1 = f1i * f2i1 * eta2all(i);
      fip = fip0 + fip1;
        
      weightjk = w(k)/M_PI;
        
      for(t = 0; t < beta1.length(); ++t){
          scores(i,t) += weightjk * Covariates(i,t) * (Y1star(i)-eta1all(i)) / pow(sigma,2) * fip;
        }
      for(t = 0; t< beta2.length(); ++t){
          scores(i,indexbeta1+t) += weightjk * Covariates(i,t) * ((0  - eta2all(i)) * fip0
                                    + (1 - eta2all(i)) * fip1) ;
        }
      for(t = 0; t < alpha1.length(); ++t){
          scores(i,indexbeta2+t) += weightjk * CovMis2(i,t) * ((1-Y2star(i)) - mta_21expit(i)) * fip1;
        }
      for(t = 0; t< alpha0.length(); ++t){
          scores(i,indexalpha1+t) += weightjk * CovMis2(i,t) * (Y2star(i) - mta_20expit(i)) * fip0;
        }
      scores(i,indexsigma-1) += weightjk * (-1/sigma + pow((Y1star(i) - eta1all(i)),2)/pow(sigma,3))* fip;
      scores(i,indexsigmag-1) += weightjk * (-1/sigma_g + uk*uk/pow(sigma_g,3)/R(i)) * fip;
      fkkt(i) +=  weightjk*fip;
      
    }
  }
   
   
  for (i = 0; i < nsam; ++i){
    for (j = 0; j < indexsigmag; ++j){
      scoresum(j) += scores(i,j) / fkkt(i);
    }
  }
  
  return(scoresum);
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
