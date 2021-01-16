# include <Rcpp.h>
using namespace Rcpp;

double fZplusYZminusj(int Yistar, int Yi, int Zminus, double muZplus, double muZminus){
  double out=0;
  
  out = R::dpois(Yistar-Yi+Zminus, exp(muZplus),0) * 
    R::dbinom(Zminus, Yi, exp(muZminus)/(1+exp(muZminus)), 0);
  
  return out;
}

// [[Rcpp::export]]

NumericVector ZI_GenerateZpZmJoint(int Yistar, int Yi, double muZplusi, double muZminusi){
  double u,P;
  int Zminusi,minindex=0,i;
  int index = 1;
  int ncandidate;
  
  // std::cout << 'P' << std::max(0,Yi-Yistar)  << ' ';
  const int nrecord = Yi - std::max(0,Yi-Yistar) + 2;
  
  NumericVector corPlus(nrecord),corMinus(nrecord),Psum(nrecord),out(2);
  
  Psum[0] = 0;
  u = R::runif(0,1);
  
  
  
  for (Zminusi = std::max(0,Yi-Yistar); Zminusi < Yi+1; ++ Zminusi){
    corMinus[index] = Zminusi;
    P = fZplusYZminusj(Yistar, Yi, Zminusi, muZplusi, muZminusi);
    Psum[index] =  Psum[index-1] + P;
    ++index;
  }
  
  // return(corMinus);
  
  ncandidate = index;

  for (i = 1; i < ncandidate; ++ i){
    if (u <= Psum[i]/Psum[ncandidate-1]) {
      minindex = i;
      break;
    }
  }


  out(0) = corMinus[minindex]; //Zminus
  out(1) = Yistar - Yi + out(0) ; // Zplus

  return out;
}


/*** R
# set.seed(2018)
# fY(3,0.5,0.5)
# result <- lapply(0:50,FUN=function(x){
#   return(ZI_GenerateBigJoint(x, 20, 20, 0.5,0.5,0.5,0.5))
# })

# ZI_GenerateJoint(6,20,0.5,3,0.5,0.5)
# ZI_GenerateZpZmJoint(2,1,-0.9273334,-3.050121)
# ZI_GenerateZpZmJoint(Ystarval[t], Yval[t], muZplus[t], muZminus[t])

# Covarplusall <- datasim$main[,c("intercept","Xplus")] #rbind(Covarplus, Covarvalplus)
# Covarminusall <- datasim$main[,c("intercept","Xminus")] #rbind(Covarminus, Covarvalminus)
# 
# muZplus <- as.matrix(Covarplusall) %*% t(t(alphaplus.candidate[,1])) 
# muZminus <- as.matrix(Covarminusall) %*% t(t(alphaminus.candidate[,1]))
# Ystarval = datasim$val$Ystar
# Yval = datasim$val$Y
# 
# t <- 8 
# ZI_GenerateZpZmJoint(Ystarval[t], Yval[t], muZplus[t], muZminus[t])

*/
