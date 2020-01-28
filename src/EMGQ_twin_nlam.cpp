# include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List EMGQ_twin_nlam(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2, double rho,  
                                    NumericVector beta1, NumericVector beta2, double gamma, NumericVector alpha1, NumericVector alpha0,
                                    double sigma, double sigma_e, double sigma_g){
  int i,j,k,l,t,s,r,m;
  int Y2;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  int ncovMis=CovMis2.ncol();
  NumericVector Y1(nsam),mta_1(nsam);
  NumericVector mta_20(nsam),mta_21(nsam),mta_20expit(nsam),mta_21expit(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam/2);
  double eta2exp,random1,random2,uk1,uk2;
  double f2i1,f2i0,f1i1,f1i0,weightjk;
  double fip_p1,fip_p2,fip0_p1,fip0_p2,fip1_p1,fip1_p2;
  double mta_11_p1,mta_10_p1,mta_11_p2,mta_10_p2;
  double Cov_meas11_p1,Cov_meas10_p1,Cov_meas11_p2,Cov_meas10_p2;
  double ui11,ui12,ui2;
  
  const int indexbeta1 = beta1.length();
  const int indexbeta2 = indexbeta1+beta2.length();
  const int indexgamma = indexbeta2+1;
  const int indexalpha1 = indexgamma+alpha1.length();
  const int indexalpha0 = indexalpha1+alpha0.length();
  const int indexsigma = indexalpha0+1;
  const int indexsigmae = indexsigma+1;
  const int indexsigmag = indexsigmae+1;
  
  const double weightPi = pow(M_PI,-2.5);
  
  NumericMatrix FisherInf_sum(indexsigmag,indexsigmag),FisherInfo_i(indexsigmag,indexsigmag);
  NumericVector scores(indexsigmag),scoresum(indexsigmag),score1_p1(indexsigmag-1),score0_p1(indexsigmag-1),score1_p2(indexsigmag-1),score0_p2(indexsigmag-1);
  double scoreg;
  
  sigma = fabs(sigma);
  sigma_e = fabs(sigma_e);
  sigma_g = fabs(sigma_g);
  
  for(i = 0; i < nsam/2; ++i)
  { 
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
  
  for(i = 0; i < nsam/2; ++i){
    for (t = 0; t < indexsigmag; ++t){
      scores(t) = 0;
      for (s = 0; s < indexsigmag; ++s){
        FisherInfo_i(t,s) = 0;
      }
    }
    
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
              
              // Person 1
              
              Y1(i) = eta1all(i) + random1;
              
              Y2 = 1;
              if (isY2in == 1){
                Cov_meas11_p1 =  Y2;
              } else {
                Cov_meas11_p1 = (2 * Y2 - 1);
              }
              mta_11_p1 = Y1(i) + gamma * Cov_meas11_p1; 
              f1i1 = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-mta_11_p1)/fabs(sigma_e), 2.0 ));
              Y2 = 0;
              if (isY2in == 1){
                Cov_meas10_p1 =  Y2;
              } else {
                Cov_meas10_p1 = (2 * Y2 - 1);
              }
              mta_10_p1= Y1(i) + gamma * Cov_meas10_p1; 
              f1i0 = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-mta_10_p1)/fabs(sigma_e), 2.0 ));
              
              f2i1 = pow(mta_21expit(i), 1-Y2star(i)) * pow(1-mta_21expit(i), Y2star(i)) ;
              f2i0 = pow(mta_20expit(i), Y2star(i)) * pow(1-mta_20expit(i), 1-Y2star(i)) ;
              
              fip0_p1 = f1i0 * f2i0 * (1-eta2all(i));
              fip1_p1 = f1i1 * f2i1 * eta2all(i);
              fip_p1 = fip0_p1 + fip1_p1;
              
              // Person 2
              
              Y1(i+nsam/2) = eta1all(i+nsam/2) + random2;
              
              Y2 = 1;
              if (isY2in == 1){
                Cov_meas11_p2 =  Y2;
              } else {
                Cov_meas11_p2 = (2 * Y2 - 1);
              }
              mta_11_p2 = Y1(i+nsam/2) + gamma * Cov_meas11_p2; 
              f1i1 = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i+nsam/2)-mta_11_p2)/fabs(sigma_e), 2.0 ));
              Y2 = 0;
              if (isY2in == 1){
                Cov_meas10_p2 =  Y2;
              } else {
                Cov_meas10_p2 = (2 * Y2 - 1);
              }
              mta_10_p2= Y1(i+nsam/2) + gamma * Cov_meas10_p2; 
              f1i0 = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i+nsam/2)-mta_10_p2)/fabs(sigma_e), 2.0 ));
              
              f2i1 = pow(mta_21expit(i+nsam/2), 1-Y2star(i+nsam/2)) * pow(1-mta_21expit(i+nsam/2), Y2star(i+nsam/2)) ;
              f2i0 = pow(mta_20expit(i+nsam/2), Y2star(i+nsam/2)) * pow(1-mta_20expit(i+nsam/2), 1-Y2star(i+nsam/2)) ;
              
              fip0_p2 = f1i0 * f2i0 * (1-eta2all(i+nsam/2));
              fip1_p2 = f1i1 * f2i1 * eta2all(i+nsam/2);
              fip_p2 = fip0_p2 + fip1_p2;
              
              for(t = 0; t < beta1.length(); ++t){
                score1_p1(t) = Covariates(i,t) * random1 / pow(sigma,2);
                score0_p1(t) = score1_p1(t);
                score1_p2(t) = Covariates(i+nsam/2,t) * random2 / pow(sigma,2);
                score0_p2(t) = score1_p2(t);
              }
              for(t = 0; t< beta2.length(); ++t){
                score1_p1(indexbeta1+t) = Covariates(i,t) * (1 - eta2all(i)) ;
                score0_p1(indexbeta1+t) = Covariates(i,t) * (0 - eta2all(i)) ;
                score1_p2(indexbeta1+t) = Covariates(i+nsam/2,t) * (1 - eta2all(i+nsam/2)) ;
                score0_p2(indexbeta1+t) = Covariates(i+nsam/2,t) * (0 - eta2all(i+nsam/2)) ;
              }
              score1_p1(indexgamma-1) =  ((Y1star(i)  - mta_11_p1)* Cov_meas11_p1)  / pow(sigma_e,2);
              score0_p1(indexgamma-1) =  ((Y1star(i)  - mta_10_p1)* Cov_meas10_p1)  / pow(sigma_e,2);
              score1_p2(indexgamma-1) =  ((Y1star(i+nsam/2)  - mta_11_p2)* Cov_meas11_p2)  / pow(sigma_e,2);
              score0_p2(indexgamma-1) =  ((Y1star(i+nsam/2)  - mta_10_p2)* Cov_meas10_p2)  / pow(sigma_e,2);
              for(t = 0; t < alpha1.length(); ++t){
                score1_p1(indexgamma+t) = CovMis2(i,t) * ((1-Y2star(i)) - mta_21expit(i)) ;
                score0_p1(indexgamma+t) = 0;
                score1_p2(indexgamma+t) = CovMis2(i+nsam/2,t) * ((1-Y2star(i+nsam/2)) - mta_21expit(i+nsam/2)) ;
                score0_p2(indexgamma+t) = 0;
              }
              for(t = 0; t< alpha0.length(); ++t){
                score1_p1(indexalpha1+t) =  0;
                score0_p1(indexalpha1+t) =  CovMis2(i,t) * (Y2star(i) - mta_20expit(i));
                score1_p2(indexalpha1+t) =  0;
                score0_p2(indexalpha1+t) =  CovMis2(i+nsam/2,t) * (Y2star(i+nsam/2) - mta_20expit(i+nsam/2));
              }
              //sigma
              score1_p1(indexsigma-1) =  -1/sigma + pow((Y1(i) - eta1all(i)),2)/pow(sigma,3);
              score0_p1(indexsigma-1) =  score1_p1(indexsigma-1);
              score1_p2(indexsigma-1) =  -1/sigma + pow((Y1(i+nsam/2) - eta1all(i+nsam/2)),2)/pow(sigma,3);
              score0_p2(indexsigma-1) =  score1_p2(indexsigma-1);
              
              //sigma_e
              score1_p1(indexsigmae-1) =  -1/sigma_e +  pow((Y1star(i) - mta_11_p1),2)/pow(sigma_e,3);
              score0_p1(indexsigmae-1) =  -1/sigma_e +  pow((Y1star(i) - mta_10_p1),2)/pow(sigma_e,3);
              score1_p2(indexsigmae-1) =  -1/sigma_e +  pow((Y1star(i+nsam/2) - mta_11_p2),2)/pow(sigma_e,3);
              score0_p2(indexsigmae-1) =  -1/sigma_e +  pow((Y1star(i+nsam/2) - mta_10_p2),2)/pow(sigma_e,3);
              
              //sigma_g
              scoreg =  -3/sigma_g + (ui11*ui11+ui12*ui12)/pow(sigma_g,3)/(1-rho) + ui2*ui2/pow(sigma_g,3)/rho;
              
              fkkt(i) +=  weightjk * fip_p1 * fip_p2;
              
              
              for(t = 0; t < indexsigmag; ++t){
                if (t < indexsigmag -1 ){
                   scores(t) += weightjk * ((score1_p1(t) * fip1_p1 + score0_p1(t) * fip0_p1) * fip_p2
                                         + (score1_p2(t) * fip1_p2 + score0_p2(t) * fip0_p2) * fip_p1);
                  
                  for (m = 0; m < indexsigmag-1; ++m){
                    FisherInfo_i(t,m) += weightjk * ((score1_p1(t)+score1_p2(t))*(score1_p1(m)+score1_p2(m)) * fip1_p1 * fip1_p2
                                                + (score1_p1(t)+score0_p2(t))*(score1_p1(m)+score0_p2(m)) * fip1_p1 * fip0_p2
                                                + (score0_p1(t)+score1_p2(t))*(score0_p1(m)+score1_p2(m)) * fip0_p1 * fip1_p2
                                                + (score0_p1(t)+score0_p2(t))*(score0_p1(m)+score0_p2(m)) * fip0_p1 * fip0_p2);
                   // FisherInfo_i(t,m) += weightjk * ((score1_p1(t) * score1_p1(m) * fip1_p1 +   score0_p1(t) * score0_p1(m) * fip0_p1)*fip_p2
                   //                    +  (score1_p2(t) * score1_p2(m) * fip1_p2 +   score0_p2(t) * score0_p2(m) * fip0_p2)*fip_p1) ;
                  }
                  FisherInfo_i(t,indexsigmag-1) += weightjk * ((score1_p1(t)+score1_p2(t))* fip1_p1 * fip1_p2
                                  + (score1_p1(t)+score0_p2(t))* fip1_p1 * fip0_p2
                                  + (score0_p1(t)+score1_p2(t)) * fip0_p1 * fip1_p2
                                  + (score0_p1(t)+score0_p2(t)) * fip0_p1 * fip0_p2) * scoreg;
                   // FisherInfo_i(t,indexsigmag-1) += weightjk * ((score1_p1(t) * fip1_p1 +   score0_p1(t) * fip0_p1)*fip_p2
                   //                             +  (score1_p2(t) * fip1_p2 +   score0_p2(t) * fip0_p2)*fip_p1) * scoreg;
                } else {
                  scores(t) += weightjk * scoreg * fip_p1 * fip_p2;
                  for (m = 0; m < indexsigmag-1; ++m){
                    FisherInfo_i(t,m) += weightjk *  ((score1_p1(m)+score1_p2(m))* fip1_p1 * fip1_p2
                                                   + (score1_p1(m)+score0_p2(m)) * fip1_p1 * fip0_p2
                                                   + (score0_p1(m)+score1_p2(m)) * fip0_p1 * fip1_p2
                                                   + (score0_p1(m)+score0_p2(m)) * fip0_p1 * fip0_p2) * scoreg;
                  }
                  FisherInfo_i(t,m) += weightjk * scoreg * scoreg * fip_p1 * fip_p2;
                }
              }
            }
            
          }
        }
      }
    }
    
    for (t = 0; t < indexsigmag; ++t){
      scoresum(t) += scores(t) / fkkt(i);
      for (m = 0; m < indexsigmag; ++m){
        FisherInf_sum(t,m) += FisherInfo_i(t,m) / fkkt(i);
      }
    }
  }
   
   
  return Rcpp::List::create(Rcpp::Named("scoresum") = scoresum,
                            Rcpp::Named("FisherInf_sum") = FisherInf_sum);
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
