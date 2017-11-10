#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double celdaCGcalcGibbsProbZ(double mCPbyS, NumericVector nCPbyTS, 
                             NumericVector nCP, int L, double alpha, double beta) {

  // Calculate for "Theta" component
  double thetaLL = log(mCPbyS + alpha);
  
  // Calculate for "Phi" component

  double b = sum(lgamma(nCPbyTS + beta));
  double d = -sum(lgamma(nCP + (L * beta)));
  
  double phiLL = b + d;
  
  double final = thetaLL + phiLL;
  return final;
}

// [[Rcpp::export]]
double celdaCGcalcGibbsProbY(NumericVector nCPbyTS, NumericVector nByTS, NumericVector nGbyTS, 
                             int nGinY, double beta, double delta, double gamma) {
  
  // Calculate for "Phi" component
  double phiLL = sum(lgamma(nCPbyTS + beta));
  
  // Calculate for "Psi" component
  double a = sum(lgamma(nGbyTS * delta));
  double d = -sum(lgamma(nByTS + (nGbyTS * delta)));

  double psiLL = a + d;
  
  // Calculate for "Eta" side
  double etaLL = log(nGinY + gamma);

  double final = phiLL + psiLL + etaLL;
  return(final);
}
