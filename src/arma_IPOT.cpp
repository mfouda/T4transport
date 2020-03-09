#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
Rcpp::List cpp_IPOT(arma::vec mu, arma::vec nu, arma::mat C, double beta, int L, int maxiter, double abstol, bool printer){
  // parameters
  int m = mu.n_elem;
  int n = nu.n_elem;
  
  // preparation
  arma::vec b(n,fill::zeros);
  for (int i=0;i<n;i++){
    b(i) = 1.0/(static_cast<double>(n));
  }
  arma::vec a(m,fill::zeros);
  
  arma::mat G = exp(-C/beta);
  arma::mat Gammaold(m,n,fill::ones);
  arma::mat Gammanew(m,n,fill::zeros);
  arma::mat Q(m,n,fill::zeros);
  double Ginc = 0.0;
  
  // main iteration
  for (int it=0;it<maxiter;it++){
    // update
    Q = G%Gammaold;
    // inner loop
    for (int l=0;l<L;l++){
      a = mu/(Q*b);
      b = nu/((Q.t())*a);
    }
    // update Gamma
    Gammanew = arma::diagmat(a)*Q*arma::diagmat(b);
    Ginc     = arma::norm(Gammaold-Gammanew,"fro");
    Gammaold = Gammanew;
    // breaking option
    if ((Ginc < abstol)&&(it>1)){
      break;
    }
    if (printer == true){
      Rcpp::Rcout << "* IPOT : iteration " << it << "/" << maxiter << " complete.." << std::endl;
    }
  }
  
  // Return
  return Rcpp::List::create(
    Rcpp::Named("plan") = Gammaold,
    Rcpp::Named("distance") = arma::accu(Gammaold%C)
  );
}