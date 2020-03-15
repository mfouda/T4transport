#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat cpp_Sinkhorn_plan(arma::vec a, arma::vec b, arma::mat costm, double lambda, int maxiter, double abstol){
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
  return(plan);
}


// [[Rcpp::export]]
Rcpp::List cpp_EPOT(arma::vec mu, arma::vec nu, arma::mat C, double beta, int L, int maxiter, double abstol, bool printer){
  // parameters
  int m = mu.n_elem;
  int n = nu.n_elem;
  
  // preparation
  arma::mat Gammaold(m,n,fill::ones);
  arma::mat Gammanew(m,n,fill::zeros);
  arma::mat Q(m,n,fill::zeros);
  double Ginc = 0.0;
  
  // main iteration
  for (int it=0;it<maxiter;it++){
    // temporary definition for C'
    Q = C - beta*arma::log(Gammaold);
    // update Gamma
    Gammanew = cpp_Sinkhorn_plan(mu,nu,Q,beta,200,abstol);
    // updating information
    Ginc     = arma::norm(Gammaold-Gammanew,"fro");
    Gammaold = Gammanew;
    // breaking option
    if ((Ginc < abstol)&&(it>1)){
      break;
    }
    if (printer == true){
      Rcpp::Rcout << "* Exact Proximal Algorithm : iteration " << it << "/" << maxiter << " complete.." << std::endl;
    }
  }
  
  // Return
  return Rcpp::List::create(
    Rcpp::Named("plan") = Gammaold,
    Rcpp::Named("distance") = arma::accu(Gammaold%C)
  );
}

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
      Rcpp::Rcout << "* Inexact Proximal Algorithm : iteration " << it << "/" << maxiter << " complete.." << std::endl;
    }
  }
  
  // Return
  return Rcpp::List::create(
    Rcpp::Named("plan") = Gammaold,
    Rcpp::Named("distance") = arma::accu(Gammaold%C)
  );
}