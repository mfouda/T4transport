#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec div_2_vecs(arma::vec x, arma::vec y){
  int m = x.n_elem;
  arma::vec output(m,fill::zeros);
  for (int i=0;i<m;i++){
    output(i) = x(i)/y(i);
  }
  return(output);
}
arma::mat compute_Gamma(arma::mat G, arma::vec x, arma::vec y){
  int m = G.n_rows;
  int n = G.n_cols;
  
  arma::mat output(m,n,fill::zeros);
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      output(i,j) = x(i)*G(i,j)*y(j);
    }
  }
  return(output);
}

// following notations for Xie (2019) "A Fast Proximal Point Method"
// [[Rcpp::export]]
arma::mat cpp_sinkhorn(arma::mat G, arma::vec wx, arma::vec wy,
                       int maxiter, double abstol){
  // parameters
  int m = G.n_rows;
  int n = G.n_cols;
  
  // variables
  arma::vec bold(n, fill::zeros);
  arma::vec bnew(n, fill::zeros);
  arma::vec aold(m, fill::zeros);
  arma::vec anew(m, fill::zeros);
  arma::mat Gammaold(m,n,fill::zeros);
  arma::mat Gammanew(m,n,fill::zeros);
  for (int i=0;i<n;i++){
    bold(i) = 1.0/(static_cast<double>(n));
  }
  
  // iteration
  double incgamma = 100000;
  int it = 0;
  arma::mat Gt = arma::trans(G);
  while (incgamma > abstol){
    // update iteration
    it += 1;
    // update a and b
    anew = div_2_vecs(wx, G*bold);
    bnew = div_2_vecs(wy, Gt*anew);
    // update Gamma & cost
    Gammanew = compute_Gamma(G, anew, bnew);
    incgamma = arma::norm(Gammaold-Gammanew,"fro");
    
    // update to old ones
    aold = anew;
    bold = bnew;
    Gammaold = Gammanew;
    if (it > maxiter){
      break;
    }
  }
  
  // return
  return(Gammaold);
}
  