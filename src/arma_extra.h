#ifndef _T4transport_ARMA_EXTRA_H
#define _T4transport_ARMA_EXTRA_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

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
// 9. extra_SubgradPlan        : compute Subgradient as well as Plan

double extra_pnorm(arma::rowvec x, arma::rowvec y, double p);
arma::mat extra_dist2mat(arma::mat X, arma::mat Y, double p);
arma::field<arma::mat> extra_dist2list(arma::mat X, arma::field<arma::mat> dlist, double p);
arma::mat extra_distmat(arma::mat &X, double p=2.0);
arma::vec extra_vec_centering(arma::vec x);
arma::vec extra_smooth_Subgradient(arma::vec a, arma::vec b, arma::mat M, double regpar, int maxiter, double abstol);
arma::vec extra_L1_normalize(arma::vec x);
bool extra_any_nan(arma::vec x);
Rcpp::List extra_SubgradPlan(arma::vec a, arma::vec b, arma::mat M, double regpar, int maxiter, double abstol);

#endif