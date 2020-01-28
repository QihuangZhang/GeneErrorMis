# include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List EMGQ_nullmeas_ci(NumericVector w, NumericVector v, NumericVector Y1star, NumericVector Y2star, NumericMatrix Covariates, int isY2in, NumericMatrix CovMis2, NumericVector R,  
                                    NumericVector beta1, NumericVector beta2, NumericVector lambda, double gamma, NumericVector alpha1, NumericVector alpha0,
                                    double sigma, double sigma_e, double sigma_g){
  int i,j,k,t,s;
  int Y2;
  const int nS= w.length();
  int nsam=Covariates.nrow();
  int ncov=Covariates.ncol();
  int ncovMis=CovMis2.ncol();
  NumericVector Y1(nsam),mta_1(nsam);
  NumericVector mta_20(nsam),mta_21(nsam),mta_20expit(nsam),mta_21expit(nsam);
  NumericVector eta1(nsam),eta2(nsam),eta1all(nsam),eta2all(nsam);
  NumericVector fkkt(nsam);
  double eta2exp,Cov_meas11,Cov_meas10,mta_11,mta_10,random1,uk;
  double fip,fip0,fip1,f2i1,f2i0,f1i1,f1i0,weightjk;
  
  const int indexbeta1 = beta1.length();
  const int indexbeta2 = indexbeta1+beta2.length();
  const int indexlambda = indexbeta2+lambda.length();
  const int indexgamma = indexlambda+1;
  const int indexalpha1 = indexgamma+alpha1.length();
  const int indexalpha0 = indexalpha1+alpha0.length();
  const int indexsigma = indexalpha0+1;
  const int indexsigmae = indexsigma+1;
  const int indexsigmag = indexsigmae+1;
  
  NumericMatrix FisherInf_sum(indexsigmag,indexsigmag),FisherInfo_i(indexsigmag,indexsigmag);
  NumericVector scores(indexsigmag),scoresum(indexsigmag),score0(indexsigmag),score1(indexsigmag);
  
  sigma = fabs(sigma);
  sigma_e = fabs(sigma_e);
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
    for (t = 0; t < indexsigmag; ++t){
      scores(t) = 0;
      for (s = 0; s < indexsigmag; ++s){
        FisherInfo_i(t,s) = 0;
      }
    }
    
    for(k = 0; k < nS; ++k){
      uk = sqrt(2*R(i)) * sigma_g * v(k);
      eta1all(i) = eta1(i) + lambda(0) * uk;
      eta2exp = eta2(i) + lambda(1) * uk;
      eta2all(i) = exp(eta2exp) / (1 + exp(eta2exp));
      if (eta2all(i) != eta2all(i)) { eta2all(i)=1 ;}
      
      for(j = 0; j < nS; ++j){
        random1 = sqrt(2) * sigma * v(j);
        Y1(i) = eta1all(i) + random1;
        
        Y2 = 1;
        if (isY2in == 1){
          Cov_meas11 =  Y2;
        } else {
          Cov_meas11 = (2 * Y2 - 1);
        }
        mta_11 = Y1(i) + gamma * Cov_meas11; 
        f1i1 = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-mta_11)/fabs(sigma_e), 2.0 ));
        Y2 = 0;
        if (isY2in == 1){
          Cov_meas10 =  Y2;
        } else {
          Cov_meas10 = (2 * Y2 - 1);
        }
        mta_10= Y1(i) + gamma * Cov_meas10; 
        f1i0 = ( 1 / (fabs(sigma_e) * sqrt(2*M_PI)) )* exp( -0.5 * pow( (Y1star(i)-mta_10)/fabs(sigma_e), 2.0 ));
        
        f2i1 = pow(mta_21expit(i), 1-Y2star(i)) * pow(1-mta_21expit(i), Y2star(i)) ;
        f2i0 = pow(mta_20expit(i), Y2star(i)) * pow(1-mta_20expit(i), 1-Y2star(i)) ;
        
        fip0 = f1i0 * f2i0 * (1-eta2all(i));
        fip1 = f1i1 * f2i1 * eta2all(i);
        fip = fip0 + fip1;
        
        for(t = 0; t < beta1.length(); ++t){
          score1(t) = Covariates(i,t) * random1 / pow(sigma,2);
          score0(t) = score1(t);
        }
        for(t = 0; t< beta2.length(); ++t){
          score1(indexbeta1+t) = Covariates(i,t) * (1 - eta2all(i)) ;
          score0(indexbeta1+t) = Covariates(i,t) * (0 - eta2all(i)) ;
        }
        score1(indexbeta2) = uk * random1/ pow(sigma,2);
        score0(indexbeta2) = score1(indexbeta2);
        score1(indexlambda-1) = uk * (1 - eta2all(i));
        score0(indexlambda-1) = uk * (0 - eta2all(i));
        score1(indexgamma-1) =  ((Y1star(i)  - mta_11)* Cov_meas11)  / pow(sigma_e,2);
        score0(indexgamma-1) =  ((Y1star(i)  - mta_10)* Cov_meas10)  / pow(sigma_e,2);
        for(t = 0; t < alpha1.length(); ++t){
          score1(indexgamma+t) = (1-Y2star(i)) - mta_21expit(i);
          score0(indexgamma+t) = 0;
        }
        for(t = 0; t< alpha0.length(); ++t){
          score1(indexalpha1+t) =  0;
          score0(indexalpha1+t) =  Y2star(i) - mta_20expit(i);
        }
        score1(indexsigma-1) =  -1/sigma + pow((Y1(i) - eta1all(i)),2)/pow(sigma,3);
        score0(indexsigma-1) =  score1(indexsigma-1);
        score1(indexsigmae-1) =  -1/sigma_e +  pow((Y1star(i) - mta_11),2)/pow(sigma_e,3);
        score0(indexsigmae-1) =  -1/sigma_e +  pow((Y1star(i) - mta_10),2)/pow(sigma_e,3);
        score1(indexsigmag-1) =  -1/sigma_g + uk*uk/pow(sigma_g,3)/R(i);
        score0(indexsigmag-1) =  score1(indexsigmag-1) ;
        
        weightjk = w(k)*w(j)/M_PI;
        fkkt(i) +=  weightjk*fip;
        for(t = 0; t < indexsigmag; ++t){
          scores(t) += weightjk * (score1(t) * fip1 + score0(t) * fip0);
          for (s = 0; s < indexsigmag; ++s){
            FisherInfo_i(t,s) += weightjk * (score1(t) * score1(s) * fip1 +
                                           score0(t) * score0(s) * fip0);
          }
        }

        }
    }
    
    for (t = 0; t < indexsigmag; ++t){
      scoresum(t) += scores(t) / fkkt(i);
      for (s = 0; s < indexsigmag; ++s){
        FisherInf_sum(t,s) += FisherInfo_i(t,s) / fkkt(i);
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
