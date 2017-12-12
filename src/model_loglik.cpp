#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

/*
*  Limited-case version of the R function. Assumes reorder=T
*/
// [[Rcpp::depends( RcppArmadillo )]]
// [[Rcpp:export]]
arma::mat rowsum(const arma::mat& countMatrix, arma::uvec group) {
 arma::mat reorderedGroupedRowSum = arma::mat(countMatrix.n_rows,
                                              countMatrix.n_cols);
 
 //TODO: Gross nested for...
 for (int groupNum=1; groupNum <= max(group); groupNum++) {
   arma::uvec groupIdx = arma::find(group == groupNum);
   reorderedGroupedRowSum[groupNum-1] = accu(countMatrix.rows(groupIdx));
 }
 
 return reorderedGroupedRowSum;
}


/*
double fastLoglikFromVariables(const arma::mat& counts, NumericVector s,
                              NumericVector z, NumericVector y,
                              int K, int L, double alpha, double beta,
                              double delta,  double gamma,
                              Function rowsum) {

   // ###################################### LAND OF C++ ########################### 
   // Calculate for "Theta" component
   NumericMatrix m  = twoFactorContingencyTable(z, s);
   int ns = m.ncol();

   double a = ns * lgamma(K * alpha);
   double b = sum(lgamma(m + alpha));
   double c = -ns * K * lgamma(alpha);
   double d = -sum(lgamma(colSums(m + alpha)));
   double thetaLL = a + b + c + d;

   //Calculate for "Phi" component
   //TODO This guy suckssss...
   arma::mat groupedRowSums = rowsum(counts, y, true).t();
   NumericVector nCPbyTS = rowsum(groupedRowSums, z);

   a = K * lgamma(L * beta);
   b = sum(lgamma(nCPbyTS + beta));
   c = -K * L * lgamma(beta);
   d = -sum(lgamma(rowSums(nCPbyTS + beta)));
   double phiLL = a + b + c + d;

    //###################################### LAND OF R ############################# 

   //## Calculate for "Psi" component
   NumericVector nByG = rowSums(counts);
   NumericVector nByTS = rowsum(nByG, y);

   NumericVector nGByTS = table(y);
   int nG = y.length();

   a = sum(lgamma(nGByTS * delta));
   b = sum(lgamma(nByG + delta));
   c = -nG * lgamma(delta);
   d = -sum(lgamma(nByTS + (nGByTS * delta)));

   double psiLL = a + b + c + d;

   //## Calculate for "Eta" side
   a = lgamma(L * gamma);
   b = sum(lgamma(nGByTS + gamma));
   c = -L * lgamma(gamma);
   d = -lgamma(sum(nGByTS + gamma));

   double etaLL = a + b + c + d;

   double result = thetaLL + phiLL + psiLL + etaLL;
   return(result);
} */
