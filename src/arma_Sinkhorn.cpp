#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// following the MATLAB code from 
// https://marcocuturi.net/SI.html
// [[Rcpp::export]]
Rcpp::List cpp_Sinkhorn_Cuturi(arma::vec a, arma::vec b, arma::mat costm, double lambda, int maxiter, double abstol){
  // prepare
  int M = costm.n_rows;
  int N = costm.n_cols;
  arma::mat K  = arma::exp(-costm/lambda);
  arma::mat Kt = arma::trans(K);

  // initialize others
  arma::mat Ktilde = arma::diagmat(1/a)*K;
  arma::vec u(M,fill::ones);
  for (int m=0;m<M;m++){
    u(m) = 1.0/(static_cast<double>(M));
  }
  arma::vec v;

  // iteration
  int iter = 0;
  double mdiff = 100000;
  while (iter < maxiter){
    iter += 1;
    u = 1 / (Ktilde * (b / (Kt*u)));
    if ((iter % 20) == 0){
      v = b/(Kt*u);
      u = 1/(Ktilde*v);
      mdiff = arma::norm((v%(Kt*u))-b);
    }
    if (mdiff < abstol){
      break;
    }
  }
  
  // let's return transport plan
  arma::mat plan  = arma::diagmat(u)*K*arma::diagmat(v);
  double distance = arma::sum(u%((K%costm)*v));
  return(Rcpp::List::create(
    Rcpp::Named("distance") = distance,
    Rcpp::Named("plan") = plan
  ));
}