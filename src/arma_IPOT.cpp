#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
Rcpp::List cpp_IPOT(arma::vec mu, arma::vec nu, arma::mat C, double beta, int L, int maxiter, double abstol){
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
  }
  
  // Return
  return Rcpp::List::create(
    Rcpp::Named("plan") = Gammaold,
    Rcpp::Named("distance") = arma::accu(Gammaold%C)
  );
}

arma::vec tmp_log(arma::vec myvec){
  int n = myvec.n_elem;
  arma::vec output(n,fill::zeros);
  for (int i=0;i<n;i++){
    output(i) = std::log(static_cast<float>(myvec(i)));
  }
  return(output);
}
arma::vec tmp_exp(arma::vec myvec){
  int n = myvec.n_elem;
  arma::vec output(n,fill::zeros);
  for (int i=0;i<n;i++){
    output(i) = std::exp(static_cast<float>(myvec(i)));
  }
  return(output);
}
arma::vec update_WB_q(arma::mat Ak, arma::mat Bk, arma::cube Hk, arma::vec lambda){
  // parameters
  int K = lambda.n_elem;
  int N = Ak.n_rows;
  
  arma::vec output(N,fill::zeros);
  for (int k=0;k<K;k++){
    output += lambda(k)*tmp_log(Ak.col(k)%(Hk.slice(k)*Bk.col(k)));
  }
  return(tmp_exp(output));
}
// [[Rcpp::export]]
Rcpp::List cpp_IPOT_WB(arma::mat C, arma::mat Pk, double beta, arma::vec lambdas, int L, int maxiter, double abstol){
  // some parameters
  int N = C.n_rows;
  int K = Pk.n_cols;
  double Ndbinv = 1.0/(static_cast<double>(N));
  
  // preliminary
  arma::mat G = exp(-C/beta);
  arma::cube Gammaold(N,N,K,fill::ones);
  arma::cube Gammanew(N,N,K,fill::ones);
  arma::cube Hcube(N,N,K,fill::zeros);
  double Gammainc = 0.0;
  
  arma::mat Amat(N,K,fill::zeros);
  arma::mat Bmat(N,K,fill::zeros);
  arma::vec q(N,fill::zeros);
  for (int n=0;n<N;n++){
    q(n) = Ndbinv;
    for (int k=0;k<K;k++){
      Amat(n,k) = Ndbinv;
      Bmat(n,k) = Ndbinv;
    }
  }
  
  // main iteration
  for (int it=0;it<maxiter;it++){
    // Line 7. Hk <- G%Hk
    for (int k=0;k<K;k++){
      Hcube.slice(k) = G%Gammaold.slice(k);
    }
    // Line 8-11. Iteration
    for (int l=0;l<L;l++){
      for (int k=0;k<K;k++){
        Amat.col(k) = q/(Hcube.slice(k)*Bmat.col(k));
        Bmat.col(k) = Pk.col(k)/(arma::trans(Hcube.slice(k))*Amat.col(k));
        q = update_WB_q(Amat, Bmat, Hcube, lambdas);
      }
    }
    // update Gamma
    for (int k=0;k<K;k++){
      Gammanew.slice(k) = arma::diagmat(Amat.col(k))*Hcube.slice(k)*arma::diagmat(Bmat.col(k));
    }
    // compute incremental change of Gamma
    Gammainc = 0.0;
    for (int k=0;k<K;k++){
      Gammainc += arma::norm(Gammaold.slice(k)-Gammanew.slice(k),"fro");
    }
    Gammaold = Gammanew;
    if ((it > 2)&&(Gammainc < abstol)){
      break;
    }
  }
  
  // now we have Gammaold in the cube form as well as q
  // Return
  return Rcpp::List::create(
    Rcpp::Named("plans") = Gammaold,
    Rcpp::Named("q") = q
  );
}
