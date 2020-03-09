#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat extra_distmat(arma::mat &X){
  // parameter
  int n = X.n_rows;
  int p = X.n_cols;
  
  // rowvecs for definition
  arma::rowvec xi(p,fill::zeros);
  arma::rowvec xj(p,fill::zeros);
  arma::rowvec xdiff(p,fill::zeros);
  
  // iterate
  double normval = 0.0;
  arma::mat output(n,n,fill::zeros);
  for (int i=0;i<(n-1);i++){
    xi = X.row(i);
    for (int j=(i+1);j<n;j++){
      xj = X.row(j);
      xdiff = xi-xj;
      normval = arma::norm(xdiff, 2);
      output(i,j) = normval;
      output(j,i) = normval;
    }
  } 
  return(output);
}