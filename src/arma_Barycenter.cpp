#include "RcppArmadillo.h"
#include "arma_extra.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// we have an issue with Numerical Precision 
// [[Rcpp::export]]
arma::vec C_Barycenter_Proximal(int N, int maxiter, double abstol,
                                  arma::field<arma::mat> listcXY,
                                  arma::field<arma::vec> marginals,
                                  arma::vec weights, double lambda,
                                  bool printer){
  // parameters
  double Ndb = static_cast<double>(N);
  int K = marginals.n_elem;
  
  // initialization
  arma::vec a_til(N,fill::zeros);
  arma::vec a_hat(N,fill::zeros);
  arma::vec alpha(N,fill::zeros);
  arma::vec a(N,fill::zeros);
  for (int n=0;n<N;n++){
    a_til(n) = 1.0/Ndb;
    a_hat(n) = 1.0/Ndb;
  }
  
  double t_0  = 2.0;
  double beta = 0.0;
  double ainc = 0.0;
  
  // iteration
  for (int it=0;it<maxiter;it++){
    beta = (static_cast<double>(it) + 1.0)/2.0;
    a = (1.0 - (1.0/beta))*a_hat + (1.0/beta)*a_til;
    alpha.fill(0.0);
    for (int k=0;k<K;k++){
      alpha += weights(k)*extra_smooth_Subgradient(a, marginals(k), listcXY(k), lambda, 1000, 0.0001);
    }
    a_til = arma::exp((-(t_0)*beta*alpha))%a_til;
    a_til = a_til/arma::accu(a_til);
    a_hat = (1.0 - (1.0)/beta)*a_hat + (1.0/beta)*a_til;
    ainc  = arma::norm(a_hat-a, 2);
    if (extra_any_nan(a_hat)){
      return(a);
    }
    if ((ainc < abstol)&&(it > 10)){
      return(a_hat);
    }
    if (printer==true){
      Rcpp::Rcout << "* BaryProximal (C++): iteration " << it+1 << "/" << maxiter  <<" completes.." << std::endl;  
    }
  }
  // main return
  return(a_hat);
}
// [[Rcpp::export]]
arma::vec C_Barycenter_Sinkhorn(int N, int maxiter, double abstol,
                                arma::field<arma::mat> listcXY,
                                arma::field<arma::vec> marginals,
                                arma::vec weights, double lambda,
                                bool printer){
  // parameters
  double Ndb = static_cast<double>(N);
  int K = marginals.n_elem;
  
  // initialization
  arma::vec a_til(N,fill::zeros);
  arma::vec a_hat(N,fill::zeros);
  arma::vec alpha(N,fill::zeros);
  arma::vec a(N,fill::zeros);
  for (int n=0;n<N;n++){
    a_til(n) = 1.0/Ndb;
    a_hat(n) = 1.0/Ndb;
  }
  
  double t_0  = 2.0;
  double beta = 0.0;
  double ainc = 0.0;
  
  // iteration
  for (int it=0;it<maxiter;it++){
    beta = (static_cast<double>(it) + 1.0)/2.0;
    a = (1.0 - (1.0/beta))*a_hat + (1.0/beta)*a_til;
    alpha.fill(0.0);
    for (int k=0;k<K;k++){
      alpha += weights(k)*extra_smooth_Subgradient(a, marginals(k), listcXY(k), lambda, 1000, 0.0001);
    }
    a_til = arma::exp((-(t_0)*beta*alpha))%a_til;
    a_til = a_til/arma::accu(a_til);
    a_hat = (1.0 - (1.0)/beta)*a_hat + (1.0/beta)*a_til;
    ainc  = arma::norm(a_hat-a, 2);
    if (extra_any_nan(a_hat)){
      return(a);
    }
    if ((ainc < abstol)&&(it > 10)){
      return(a_hat);
    }
    if (printer==true){
      Rcpp::Rcout << "* BarySinkhorn (C++): iteration " << it+1 << "/" << maxiter  <<" completes.." << std::endl;  
    }
  }
  // main return
  return(a_hat);
}