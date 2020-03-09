#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List cpp_Sinkhorn_Bregman(arma::vec mu, arma::vec nu, arma::mat costm, double lambda, int maxiter, double abstol, bool printer){
  // prepare
  int M = costm.n_rows;
  int N = costm.n_cols;
  arma::mat G  = arma::exp(-costm/lambda);
  
  // initialize others
  arma::vec b(N,fill::zeros);
  for (int n=0;n<N;n++){
    b(n) = 1.0/(static_cast<double>(N));
  }
  arma::vec a(M,fill::zeros);
  arma::mat Gammaold = arma::diagmat(a)*G*arma::diagmat(b);
  arma::mat Gammanew(M,N,fill::zeros);
  double    Gammainc = 0.0;
  
  // iteration
  int iter = 0;
  double mdiff = 10000;
  if (maxiter < 200){
    maxiter = 200;
  }
  while (iter < maxiter){
    iter += 1;
    a = mu / (G*b);
    b = nu / (G.t()*a);
    if (printer == true){
      Rcpp::Rcout << "* Sinkhorn : iteration " << iter << "/" << maxiter << " complete.." << std::endl;
    }
    if ((iter % 20) == 0){
      Gammanew = arma::diagmat(a)*G*arma::diagmat(b);
      Gammainc = arma::norm(Gammanew-Gammaold,"fro");
      Gammaold = Gammanew;
      if (Gammainc < abstol){
        break;
      }
    }
  }
  
  // compute and return
  return(Rcpp::List::create(
      Rcpp::Named("distance") = arma::accu(Gammaold%costm),
      Rcpp::Named("plan") = Gammaold
  ));
}



// following the MATLAB code from 
// https://marcocuturi.net/SI.html
// [[Rcpp::export]]
Rcpp::List cpp_Sinkhorn_Cuturi(arma::vec a, arma::vec b, arma::mat costm, double lambda, int maxiter, double abstol, bool printer){
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
    if (printer == true){
      Rcpp::Rcout << "* Sinkhorn : iteration " << iter << "/" << maxiter << " complete.." << std::endl;
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