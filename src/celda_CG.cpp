#include <RcppArmadillo.h>
using namespace Rcpp;


/*
*  Limited-case version of the R function. Assumes reorder=T
*/
// [[Rcpp::depends( RcppArmadillo )]]
// [[Rcpp:export]]
arma::mat rowsum(const arma::mat& countMatrix, arma::vec group) {
 arma::mat reorderedGroupedRowSum = arma::mat(countMatrix.n_rows,
                                              countMatrix.n_cols);
 for (int groupNum=1; groupNum <= max(group); groupNum++) {
   arma::uvec groupIdx = arma::find(group == groupNum);
   reorderedGroupedRowSum[groupNum-1] = accu(countMatrix.rows(groupIdx));
 }
 return reorderedGroupedRowSum;
}


NumericVector sampleLL(NumericVector llProbs) {
  NumericVector probsSub = exp(llProbs - max(llProbs));
  NumericVector probsNorm = probsSub / sum(probsSub);
  NumericVector probsSelect = sample(NumericVector(1, llProbs.length()), 1, false, probsNorm);
  return probsSelect;
}


// [[Rcpp::depends( RcppArmadillo )]]
// [[Rcpp:export]]
double celdaCGgibbsSampling(int iter, const arma::mat& counts, NumericVector y, NumericVector z,
                            double alpha, double beta,  // Dirichlet parameters
                            int K, int L, 
                            int s, // Sample Labels
                            NumericMatrix mCPbyS,  //something Cell Populations by State? cols?
                            NumericMatrix nCPbyS,  // something Cell Populations by State? rows?
                            NumericVector nCP, // number of Cell Populations
                            NumericMatrix nCPbyTS, //Cell Populations by Trscrpt. State?
                            NumericMatrix nTSbyC   // Transcriptional States by Cells
                              ) {
  
  //Start by zero-indexing any variable that's used for vector or matrix accession:
  y = y - 1;
  z = z - 1; 
  
  // Begin process of Gibbs sampling for each cell
  nTSbyC = rowsum(counts, y);
  NumericVector ix = sample(y, y.length());
  for (int i: y) { //350
    // Subtract current cell counts from matrices
    mCPbyS(z[i], s[i]) = mCPbyS(z[i], s[i]) - 1; //353
    nCPbyTS(z[i], _) = nCPbyTS(z[i], _) - nTSbyC(_, i);
    nCP[z[i]] = nCP[z[i]] - sum(counts(_, i));
    
    //Calculate probabilites for each state
    NumericVector probs(K);
    for (int j: range(0, K-1)) { //359
      NumericMatrix nCPbyTStemp(clone(nCPbyTS));
      nCPbyTStemp(j, _) = nCPbyTStemp(j, _) + nTSbyC(_, i);
      NumericVector nCPtemp(clone(nCP));
      nCPtemp[j] = nCPtemp[j] + sum(counts(_, i));
      probs[j] = celdaCGcalcGibbsProbZ(mCPbyS(j, s[i]),
                                       nCPbyTStemp, nCPtemp,
                                       L, alpha, beta); //365
    }
    
    // Sample the next state and add the counts back in
    NumericVector previousZ = z;
    z[i] = sampleLL(probs); //370
  }
  
  NumericVector endZ = z + 1;
  NumericVector endY = y + 1;
  return 666.0;
}


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
