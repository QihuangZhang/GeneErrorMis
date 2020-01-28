# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector IL_score_measonly(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2, NumericVector R,  
                                    NumericVector beta1, NumericVector beta2, double gamma,
                                    double sigma, double sigma_e, double sigma_g){
  int i,j,k,t;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  NumericVector Y1(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam);
  double eta2exp,Cov_meas11,mta_1,random1,uk;
  double fip,f2i,f1i,weightjk;
  
  const int indexbeta1 = beta1.length();
  const int indexbeta2 = indexbeta1+beta2.length();
  const int indexgamma = indexbeta2+1;
  const int indexsigma = indexgamma+1;
  const int indexsigmae = indexsigma+1;
  const int indexsigmag = indexsigmae+1;
  
  NumericMatrix scores(nsam/2,indexsigmag);
  NumericVector scoresum(indexsigmag);
  
  sigma = fabs(sigma);
  sigma_e = fabs(sigma_e);
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
  
  for(i = 0; i < nsam; ++i){
    for(k = 0; k < nS; ++k){
      uk = sqrt(2*R(i)) * sigma_g * v(k);
      eta1all(i) = eta1(i) +  uk;
      eta2exp = eta2(i) +  uk;
      eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
      if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}
      
      for(j = 0; j < nS; ++j){
        random1 = sqrt(2) * sigma * v(j);
        Y1(i) = eta1all(i) + random1;
        
        if (Y2star(i) == 1) {
          f2i = eta2all(i);
        } else {
          f2i = 1 - eta2all(i);
        }
        
         
        if (isY2in == 1){
          Cov_meas11 =  Y2star(i);
        } else {
          Cov_meas11 = (2 * Y2star(i) - 1);
        }
        mta_1 = Y1(i) + gamma * Cov_meas11; 
        f1i = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) ) * exp( -0.5 * pow( (Y1star(i)-mta_1)/fabs(sigma_e), 2.0 ));
        
        fip = f1i * f2i;
        
        weightjk = w(k)*w(j)/M_PI;
        
        for(t = 0; t < beta1.length(); ++t){
          scores(i,t) += weightjk * fip * Covariates(i,t) * random1 / pow(sigma,2);
        }
        for(t = 0; t< beta2.length(); ++t){
          scores(i,indexbeta1+t) += weightjk * fip * Covariates(i,t) * (Y2star(i)   - eta2all(i));
        }
        scores(i,indexgamma-1) += weightjk * fip * (Y1star(i)  - mta_1) * Cov_meas11  / pow(sigma_e,2);
        scores(i,indexsigma-1) += weightjk * fip * (-1/sigma + pow((Y1(i) - eta1all(i)),2)/pow(sigma,3));
        scores(i,indexsigmae-1) += weightjk * fip * (-1/sigma_e  + pow((Y1star(i) - mta_1),2)/pow(sigma_e,3));
        scores(i,indexsigmag-1) += weightjk * fip * (-1/sigma_g + uk*uk/pow(sigma_g,3)/R(i));
        fkkt(i) +=  weightjk*fip;
      }
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
