# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix IL_infomat_twin_misconly(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2,  double rho, 
                                    NumericVector beta1, NumericVector beta2, NumericVector alpha1, NumericVector alpha0,
                                    double sigma, double sigma_g){
  int i,j,k,t,l,s;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  int ncovMis=CovMis2.ncol();
  NumericVector Y1(nsam),mta_1(nsam);
  NumericVector mta_20(nsam),mta_21(nsam),mta_20expit(nsam),mta_21expit(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam/2);
  double eta2exp,uk1,uk2;
  double fip_1p,fip_2p,fip0_1p,fip0_2p,fip1_1p,fip1_2p,f2i1,f2i0,f1i,weightjk;
  
  
  double ui11,ui12,ui2;
  
  const int indexbeta1 = beta1.length();
  const int indexbeta2 = indexbeta1+beta2.length();
  const int indexalpha1 = indexbeta2+alpha1.length();
  const int indexalpha0 = indexalpha1+alpha0.length();
  const int indexsigma = indexalpha0+1;
  const int indexsigmag = indexsigma+1;
  
  
  const double weightPi = pow(M_PI,-0.5*3);
  
  NumericMatrix scores(nsam/2,indexsigmag),scoreswtd(nsam/2,indexsigmag),infoMat(indexsigmag,indexsigmag);

  
  sigma = fabs(sigma);
  sigma_g = fabs(sigma_g);
  
  for(i = 0; i < nsam/2; ++i){ 
    fkkt(i)=0;
  }
  
  for(i = 0; i < nsam; ++i)
  { eta1(i) = 0;
    eta2(i) = 0;
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
  
  for (j = 0; j < indexsigmag; ++j){
    for (k = 0; k < indexsigmag; ++k){
      infoMat(j,k) = 0;
    }
  }
  
  for(i = 0; i < nsam/2; ++i){
    for(k = 0; k < nS; ++k){
      for(l = 0; l < nS; ++l){
        for(s = 0; s < nS; ++s){
          weightjk = w(k) * w(l) * w(s) * weightPi;
          ui11 = sqrt(2*(1-rho)) * sigma_g * v(k);
          ui12 = sqrt(2*(1-rho)) * sigma_g * v(l);
          ui2 = sqrt(2*(rho)) * sigma_g * v(s);
          
          //person 1
          uk1 = ui11 + ui2;
          eta1all(i) = eta1(i) +  uk1;
          eta2exp = eta2(i) +  uk1;
          eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
          if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}
          
          f1i = ( 1 / (fabs(sigma) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-eta1all(i))/fabs(sigma), 2.0 ));
          
          f2i1 = pow(mta_21expit(i), 1-Y2star(i)) * pow(1-mta_21expit(i), Y2star(i)) ;
          f2i0 = pow(mta_20expit(i), Y2star(i)) * pow(1-mta_20expit(i), 1-Y2star(i)) ;
          
          fip0_1p = f1i * f2i0 * (1-eta2all(i));
          fip1_1p = f1i * f2i1 * eta2all(i);
          fip_1p = fip0_1p + fip1_1p;
          
          // Person 2
          uk2 = ui12 + ui2;
          eta1all(i + nsam/2) = eta1(i + nsam/2) +  uk2;
          eta2exp = eta2(i + nsam/2) +  uk2;
          eta2all(i + nsam/2) = exp(eta2exp) / (1 + exp(eta2exp));
          if (eta2all(i + nsam/2) != eta2all(i + nsam/2)) { eta2all(i + nsam/2)=1 ;}
          
          f1i = ( 1 / (fabs(sigma) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i + nsam/2)-eta1all(i + nsam/2))/fabs(sigma), 2.0 ));
          
          f2i1 = pow(mta_21expit(i + nsam/2), 1-Y2star(i + nsam/2)) * pow(1-mta_21expit(i + nsam/2), Y2star(i + nsam/2)) ;
          f2i0 = pow(mta_20expit(i + nsam/2), Y2star(i + nsam/2)) * pow(1-mta_20expit(i + nsam/2), 1-Y2star(i + nsam/2)) ;
          
          fip0_2p = f1i * f2i0 * (1-eta2all(i + nsam/2));
          fip1_2p = f1i * f2i1 * eta2all(i + nsam/2);
          fip_2p = fip0_2p + fip1_2p;
          
          
          for(t = 0; t < beta1.length(); ++t){
            scores(i,t) += weightjk * (Covariates(i,t) * (Y1star(i)-eta1all(i)) +
              Covariates(i + nsam/2,t) * (Y1star(i + nsam/2)-eta1all(i + nsam/2))) / pow(sigma,2) * fip_1p * fip_2p;}
          for(t = 0; t< beta2.length(); ++t){
            scores(i,indexbeta1+t) += weightjk * (Covariates(i,t) * ((0  - eta2all(i)) * fip0_1p + (1 - eta2all(i)) * fip1_1p) * fip_2p
                                                    + Covariates(i + nsam/2,t) * ((0  - eta2all(i + nsam/2)) * fip0_2p + (1 - eta2all(i + nsam/2)) * fip1_2p) * fip_1p) ;
          }
          for(t = 0; t < alpha1.length(); ++t){
            scores(i,indexbeta2+t) += weightjk * (CovMis2(i,t)  * ((1-Y2star(i)) - mta_21expit(i)) * fip1_1p * fip_2p
                                                    + CovMis2(i + nsam/2,t)  * ((1-Y2star(i + nsam/2)) - mta_21expit(i + nsam/2)) * fip1_2p * fip_1p);
          }
          for(t = 0; t< alpha0.length(); ++t){
            scores(i,indexalpha1+t) += weightjk * (CovMis2(i,t)  * (Y2star(i) - mta_20expit(i)) * fip0_1p * fip_2p
                                                     + CovMis2(i + nsam/2,t)  * (Y2star(i + nsam/2) - mta_20expit(i + nsam/2)) * fip0_2p * fip_1p);
          }
          
          scores(i,indexsigma-1) += weightjk * (-2/sigma + pow((Y1star(i) - eta1all(i)),2)/pow(sigma,3)
                                                  +  pow((Y1star(i + nsam/2) - eta1all(i + nsam/2)),2)/pow(sigma,3))* fip_1p * fip_2p;
          scores(i,indexsigmag-1) += weightjk * (-3/sigma_g + (ui11*ui11+ui12*ui12)/pow(sigma_g,3)/(1-rho) + ui2*ui2/pow(sigma_g,3)/rho) * fip_1p * fip_2p;
          fkkt(i) +=  weightjk * fip_1p * fip_2p;
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
