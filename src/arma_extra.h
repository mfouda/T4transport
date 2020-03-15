#ifndef _T4transport_ARMA_EXTRA_H
#define _T4transport_ARMA_EXTRA_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double extra_pnorm(arma::rowvec x, arma::rowvec y, double p);
arma::mat extra_dist2mat(arma::mat X, arma::mat Y, double p);
arma::field<arma::mat> extra_dist2list(arma::mat X, arma::field<arma::mat> dlist, double p);
arma::mat extra_distmat(arma::mat &X, double p=2.0);
arma::vec extra_vec_centering(arma::vec x);
arma::vec extra_smooth_Subgradient(arma::vec a, arma::vec b, arma::mat M, double regpar, int maxiter, double abstol);
arma::vec extra_L1_normalize(arma::vec x);
bool extra_any_nan(arma::vec x);

#endif