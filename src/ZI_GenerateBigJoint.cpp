# include <Rcpp.h>
using namespace Rcpp;

double fZplusYZminusU1U2(int Zplus, int Yi, int Zminus, int Ui1, int Ui2,
                     double  mu1, double mu2, double muZplus, double muZminus){
  double out=0;
  
  out = R::dpois(Zplus, exp(muZplus),0) *
    // R::dbinom(Zminus, Yi, R::pnorm(muZminus,0,1,1,0), 0) * 
    R::dbinom(Zminus, Yi, exp(muZminus)/(1+exp(muZminus)), 0) *
    R::dpois(Ui1,exp(mu1),0) *
    R::dpois(Ui2,exp(mu2),0);
  
  return out;
}

// [[Rcpp::export]]

NumericVector ZI_GenerateBigJoint(int Yistar, int Ui1bound, int Ui2bound, double mu1i, double mu2i, double muZplusi, double muZminusi){
  double u,P;
  int Ui1,Ui2,Zminusi,minindex=0,Yi,i;
  int index = 1;
  const int nrecord = (Yistar+1)*(Ui1bound+1)*(Ui2bound+1)+1;
  int ncandidate;

  NumericVector corPlus(nrecord),corMinus(nrecord),corU1(nrecord),corU2(nrecord),Psum(nrecord),out(5);

  Psum[0] = 0;
  u = R::runif(0,1);
  
  for (Ui1 = 0; Ui1 < Ui1bound+1; ++ Ui1){
    for (Ui2 = 0; Ui2 < Ui2bound+1; ++ Ui2){
      if (Ui1==0) {
        Yi = 0;
      } else {
        Yi = Ui2;
      }
      // std::cout << 'P' << std::max(0,Yi-Yistar)  << ' ';
      for (Zminusi = std::max(0,Yi-Yistar); Zminusi < Yi+1; ++ Zminusi){
        corMinus[index] = Zminusi;
        corU1[index] = Ui1;
        corU2[index] = Ui2;
        P = fZplusYZminusU1U2(Yistar+Zminusi-Yi, Yi, Zminusi, Ui1, Ui2,
                          mu1i, mu2i, muZplusi, muZminusi);
        Psum[index] =  Psum[index-1] + P;
        ++index;
      }
    }
  }

  ncandidate = index;
  
  for (i = 1; i < ncandidate; ++ i){
    if (u <= Psum[i]/Psum[ncandidate-1]) {
      minindex = i;
      break;
    }
  }

  out(3) = corU1[minindex];  //U1
  out(4) = corU2[minindex];  //U2
  out(1) = corMinus[minindex]; //Zminus
  if (out(3)==0) { 
    out(2) = 0;            //Y
  } else {
    out(2) = out(4);
  }
  out(0) = Yistar + out(1) - out(2); // Zplus

  return out;
}


/*** R
# set.seed(2018)
# fY(3,0.5,0.5)
# result <- lapply(0:50,FUN=function(x){
#   return(ZI_GenerateBigJoint(x, 20, 20, 0.5,0.5,0.5,0.5))
# })

# ZI_GenerateJoint(6,20,0.5,3,0.5,0.5)
# GenerateJoint(5,15,0.5,-0.7,0.5,0.3) 
*/
