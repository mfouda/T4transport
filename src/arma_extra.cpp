#include "RcppArmadillo.h"
#include "arma_extra.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// Functions
// 1. extra_pnorm              : compute |x-y|_p
// 2. extra_dist2mat           : compute d(x,y)      from {X,Y}
// 3. extra_dist2list          : compute dxy's       from {X, {Y}}
// 4. extra_distmat            : compute d(x_i, x_j) from {X}
// 5. extra_vec_centering      : take the centering of a vector
// 6. extra_smooth_Subgradient : compute subgradient
// 7. extra_L1_normalize       : do the L1 normalization of a vector  
// 8. extra_any_nan            : BOOL if a vector has any NaN/NA value


// 1. extra_pnorm : compute |x-y|_p
// [[Rcpp::export]]
double extra_pnorm(arma::rowvec x, arma::rowvec y, double p){
  // parameter
  int n = x.n_elem;
  double output = 0.0;
  
  for (int i=0;i<n;i++){
    output += std::pow(std::abs(x(i)-y(i)), p);
  }
  return(std::pow(output,1.0/p));
}

// 2. extra_dist2mat : compute d(x,y) from {X,Y}
// [[Rcpp::export]]
arma::mat extra_dist2mat(arma::mat X, arma::mat Y, double p){
  // parameters
  int m = X.n_rows;
  int n = Y.n_rows;
  int d = X.n_cols;
  
  // iteration
  arma::mat output(m,n,fill::zeros);
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      output(i,j) = extra_pnorm(X.row(i), Y.row(j), p);
    }
  }
  
  // return
  return(output);
}

// 3. extra_dist2list : compute dxy's from {X, {Y}}
// [[Rcpp::export]]
arma::field<arma::mat> extra_dist2list(arma::mat X, arma::field<arma::mat> dlist, double p){
  // parameter
  int K = dlist.n_elem;
  arma::field<arma::mat> outdist(K);
  
  for (int k=0;k<K;k++){
    outdist(k) = extra_dist2mat(X, dlist(k), p);
  }
  // return
  return(outdist);
}

// 4. extra_distmat : compute d(x_i, x_j) from {X}
// [[Rcpp::export]]
arma::mat extra_distmat(arma::mat &X, double p){
  // parameter
  int n = X.n_rows;
  int d = X.n_cols;
  
  // rowvecs for definition
  arma::rowvec xi(d,fill::zeros);
  arma::rowvec xj(d,fill::zeros);

  // iterate
  double normval = 0.0;
  arma::mat output(n,n,fill::zeros);
  for (int i=0;i<(n-1);i++){
    xi = X.row(i);
    for (int j=(i+1);j<n;j++){
      xj = X.row(j);
      normval = extra_pnorm(xi, xj, p);
      output(i,j) = normval;
      output(j,i) = normval;
    }
  } 
  return(output);
}

// 5. extra_vec_centering : take the centering of a vector
// [[Rcpp::export]]
arma::vec extra_vec_centering(arma::vec x){
  return(x-arma::mean(x));
}


// 6. extra_smooth_Subgradient : compute subgradient



// [[Rcpp::export]]
arma::vec extra_smooth_Subgradient(arma::vec a, arma::vec b, arma::mat M, double regpar, int maxiter, double abstol) {
  // prepare
  double lambda = 1.0/regpar;           // corresponding to my convention 
  arma::mat K   = arma::exp(-lambda*M);
  arma::mat Ktil = arma::diagmat(1/a)*K;
  
  // initialize
  int N       = a.n_elem;
  double Ndb  = static_cast<double>(N);
  double uinc = 0.0;
  arma::vec uold(N,fill::zeros);
  arma::vec unew(N,fill::zeros);
  uold = uold.ones()/Ndb;
  unew = unew.ones()/Ndb;
  
  
  
  // main iteration
  for (int it=0; it<maxiter; it++){
    // Sinkhorn update
    uold = 1.0/(Ktil*(b/(K.t()*uold)));
    // if (extra_any_nan(unew)){
    //   break;
    // }
    // uold = unew;
    // while u changes.. do 
    if (it%5 == 0){
      unew = 1.0/(Ktil*(b/(K.t()*uold)));
      uinc = arma::norm(uold-unew,2);
      if (uinc < abstol){
        uold = unew;
        break;
      } else {
        uold = unew;
      }
    }
  }
  
  // compute smoothed dual solution
  arma::vec alpha  = (1.0/lambda)*arma::log(uold) - arma::accu(arma::log(uold))/(lambda*Ndb);
  return(alpha);
}

// 7. extra_L1_normalize : do the L1 normalization of a vector
// [[Rcpp::export]]
arma::vec extra_L1_normalize(arma::vec x){
  double xsum = arma::accu(x);
  return(x/xsum);
}

// 8. extra_any_nan : BOOL if a vector has any NaN/NA value
// [[Rcpp::export]]
bool extra_any_nan(arma::vec x){
  int N = x.n_elem;
  for (int n=0;n<N;n++){
    if (!arma::is_finite(x(n))){
      return(true);
    }
  }
  return(false);
}