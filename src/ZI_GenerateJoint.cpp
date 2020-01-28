# include <Rcpp.h>
using namespace Rcpp;


double fY(int Yi, double phi, double mu2) {
  double out=0;
  
  if (Yi==0){
    out = 1-phi + phi * R::dpois(0,exp(mu2),0);
  } else {
    out = phi * R::dpois(Yi,exp(mu2),0);
  }
  
  return out;
}

double fZplusYZminus(int Zplus, int Yi, int Zminus,
                     double  phi, double mu2, double muZplus, double muZminus){
  double out=0;
  
  out = fY(Yi, phi, mu2) *  R::dpois(Zplus, exp(muZplus),0) *
    R::dbinom(Zminus, Yi, R::pnorm(muZminus,0,1,1,0), 0);
  
  return out;
}

// [[Rcpp::export]]

NumericVector ZI_GenerateJoint(int Yistar, int minusbound, double phii, double mu2i, double muZplusi, double muZminusi){
  double u,P;
  int Zplusi,Zminusi,minindex=0,i;
  int index = 1;
  const int nrecord = (Yistar+1)*(minusbound+1)+1;

  NumericVector corPlus(nrecord),corMinus(nrecord),Psum(nrecord),out(3);

  Psum[0] = 0;
  u = R::runif(0,1);

  for (Zplusi = 0; Zplusi < Yistar+1; ++ Zplusi){
    for (Zminusi =0; Zminusi < minusbound + 1; ++ Zminusi) {
      corPlus[index] = Zplusi;
      corMinus[index] = Zminusi;
      P = fZplusYZminus(Zplusi, Yistar-Zplusi+Zminusi, Zminusi, phii, mu2i, muZplusi, muZminusi);
      Psum[index] =  Psum[index-1] + P;
      ++index;
    }
  }

  for (i = 1; i <= nrecord; ++ i){
    if (u <= Psum[i]/Psum[nrecord-1]) {
      minindex = i;
      break;
    }
  }

  out[0] = corPlus[minindex];  // Zplus
  out[1] = corMinus[minindex];  // Zminus
  out[2] = Yistar - out(0) + out(1);  // Y

  return out;
}


/*** R
# set.seed(2018)
# fY(3,0.5,0.5)
# fZplusYZminus(3, 0, 0, 0.5,0.5,0.5,0.5)
# ZI_GenerateJoint(6,20,0.5,3,0.5,0.5)
# GenerateJoint(5,15,0.5,-0.7,0.5,0.3) 
*/
