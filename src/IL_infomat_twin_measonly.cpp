# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix IL_infomat_twin_measonly(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2, double rho, 
                                    NumericVector beta1, NumericVector beta2, double gamma, 
                                    double sigma, double sigma_e, double sigma_g){
  int i,j,k,t,l,s,r;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  NumericVector Y1(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam/2);
  double eta2exp,Cov_meas11_1p,Cov_meas11_2p,mta_11_1p,mta_11_2p,random1,random2,uk1,uk2;
  double fip_p1,fip_p2,f2i,f1i,weightjk;
  double ui11,ui12,ui2;
  
  const int indexbeta1 = beta1.length();
  const int indexbeta2 = indexbeta1+beta2.length();
  const int indexgamma = indexbeta2+1;
  const int indexsigma = indexgamma+1;
  const int indexsigmae = indexsigma+1;
  const int indexsigmag = indexsigmae+1;
  
  
  const double weightPi = pow(M_PI,-2.5);
  
  
  NumericMatrix scores(nsam/2,indexsigmag),scoreswtd(nsam/2,indexsigmag),infoMat(indexsigmag,indexsigmag);

  
  sigma = fabs(sigma);
  sigma_e = fabs(sigma_e);
  sigma_g = fabs(sigma_g);
  
  for(i = 0; i < nsam/2; ++i){ 
    fkkt(i)=0;}
  
  for(i = 0; i < nsam; ++i)
  { eta1(i) = 0;
    eta2(i) = 0;
    for (j = 0; j < ncov; ++j){
      eta1(i) += beta1(j) * Covariates(i,j);
      eta2(i) += beta2(j) * Covariates(i,j);
    }
  }
  
  for (j = 0; j < indexsigmag; ++j){
    for (k = 0; k < indexsigmag; ++k){
      infoMat(j,k) = 0;
    }
  }
  
  for(i = 0; i < nsam/2; ++i){
    for(k = 0; k < nS; ++k){
      for(l = 0; l < nS; ++l){
        for(s = 0; s < nS; ++s){
          ui11 = sqrt(2*(1-rho)) * sigma_g * v(k);
          ui12 = sqrt(2*(1-rho)) * sigma_g * v(l);
          ui2 = sqrt(2*(rho)) * sigma_g * v(s);
          
          //person 1
          uk1 = ui11 + ui2;
          eta1all(i) = eta1(i) +  uk1;
          eta2exp = eta2(i) +  uk1;
          eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
          if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}
          
          // Person 2
          uk2 = ui12 + ui2;
          eta1all(i + nsam/2) = eta1(i + nsam/2) +  uk2;
          eta2exp = eta2(i + nsam/2) +  uk2;
          eta2all(i + nsam/2) = exp(eta2exp) / (1 + exp(eta2exp));
          if (eta2all(i + nsam/2) != eta2all(i + nsam/2)) { eta2all(i + nsam/2)=1 ;}
          
          for(j = 0; j < nS; ++j){
            for(r = 0; r < nS; ++r){
              weightjk = w(k)*w(l)*w(s)*w(j)*w(r)*weightPi;
              random1 = sqrt(2) * sigma * v(j);
              random2 = sqrt(2) * sigma * v(r);
              
              //Person 1:
              Y1(i) = eta1all(i) + random1;
              
              
              if (Y2star(i) == 1) {
                f2i = eta2all(i);
              } else {
                f2i = 1 - eta2all(i);
              }
              
              
              if (isY2in == 1){
                Cov_meas11_1p =  Y2star(i);
              } else {
                Cov_meas11_1p = (2 * Y2star(i) - 1);
              }
              mta_11_1p = Y1(i) + gamma * Cov_meas11_1p; 
              f1i = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-mta_11_1p)/fabs(sigma_e), 2.0 ));
              
              fip_p1 = f1i * f2i;
              
              // Person 2:
              Y1(i+nsam/2) = eta1all(i+nsam/2) + random2;
              
              if (Y2star(i+nsam/2) == 1) {
                f2i = eta2all(i+nsam/2);
              } else {
                f2i = 1 - eta2all(i+nsam/2);
              }
              
              
              if (isY2in == 1){
                Cov_meas11_2p =  Y2star(i+nsam/2);
              } else {
                Cov_meas11_2p = (2 * Y2star(i+nsam/2) - 1);
              }
              mta_11_2p = Y1(i+nsam/2) + gamma * Cov_meas11_2p; 
              f1i = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i+nsam/2)-mta_11_2p)/fabs(sigma_e), 2.0 ));
              
              fip_p2 = f1i * f2i;
              
              for(t = 0; t < beta1.length(); ++t){
                scores(i,t) += weightjk * fip_p1 * fip_p2 * (Covariates(i,t) * random1 + Covariates(i+nsam/2,t) * random2)  / pow(sigma,2);
              }
              for(t = 0; t< beta2.length(); ++t){
                scores(i,indexbeta1+t) += weightjk * fip_p1 * fip_p2 * (Covariates(i,t) * (Y2star(i)   - eta2all(i)) 
                                                                          + Covariates(i+nsam/2,t) * (Y2star(i+nsam/2)   - eta2all(i+nsam/2)));
              }
              scores(i,indexgamma-1) += weightjk * fip_p1 * fip_p2 * ((Y1star(i)  - mta_11_1p) * Cov_meas11_1p  / pow(sigma_e,2)
                                                                        +  (Y1star(i+nsam/2)  - mta_11_2p) * Cov_meas11_2p  / pow(sigma_e,2));
              scores(i,indexsigma-1) += weightjk * fip_p1 * fip_p2 * (-2/sigma + pow((Y1(i) - eta1all(i)),2)/pow(sigma,3)
                                                                        + pow((Y1(i+nsam/2) - eta1all(i+nsam/2)),2)/pow(sigma,3));
              scores(i,indexsigmae-1) += weightjk * fip_p1 * fip_p2 * (-2/sigma_e  + pow((Y1star(i) - mta_11_1p),2)/pow(sigma_e,3)
                                                                         + pow((Y1star(i+nsam/2) - mta_11_2p),2)/pow(sigma_e,3));
              scores(i,indexsigmag-1) += weightjk * fip_p1 * fip_p2 * (-3/sigma_g + (ui11*ui11+ui12*ui12)/pow(sigma_g,3)/(1-rho) + ui2*ui2/pow(sigma_g,3)/rho);
              fkkt(i) +=  weightjk * fip_p1 * fip_p2;
            }
            
          }
        }
      }
    }
  }
   
 
   
  for (i = 0; i < nsam/2; ++i){
    for (j = 0; j < indexsigmag; ++j){
      scoreswtd(i,j) = scores(i,j) / fkkt(i);
    }
    for (j = 0; j < indexsigmag; ++j){
      for (k = 0; k < indexsigmag; ++k){
        infoMat(j,k) += scoreswtd(i,j)*scoreswtd(i,k);
      }
    }
  }
  
  return(infoMat);
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
