#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec cpp_Sinkhorn_Subgradient(arma::vec a, arma::vec b, arma::mat M, double regeps, int maxIter = 1000, double tolerance = 0.0001) {
  
  // Transforming Input, i.e. calculating the kernel
  
  double lambda = 1.0/(static_cast<double>(regeps));
  arma::mat K = exp(-lambda*M);
  arma::mat ainvK = diagmat(1/a) * K;
  
  
  //Initialize u, v and the vector next for the stopping criterion
  arma::vec u(a.n_rows);
  u = u.ones()/a.n_rows;
  arma::vec next;
  
  for(int i=0;i<maxIter;i++){
    
    u = 1 / (ainvK * (b / (K.t()*u)));  //Sinkhorn`s update
    
    //The stopping criterion is checked every 20th step, i.e. if u and the next update of u (called next) differ in 2-norm only by a tolerance
    if(i % 20 == 0){
      next = 1 / (ainvK * (b / (K.t()*u)));
      double Criterion = norm(abs(next-u));
      if(Criterion < tolerance){
        u = next;
        break;
      }
      else{
        u = next;
      }
    }
    
  }
  
  //Computation of the dual solution (c.f. Proposition 2, Fast Computation of Wasserstein Barycenters)
  
  arma::vec alpha = (1 / lambda) * log(u) - accu(log(u))/(lambda*(u.n_rows));
  
  return(alpha);
}

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