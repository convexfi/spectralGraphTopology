#ifndef LINALG_H
#define LINALG_H
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

arma::vec eigval_sym(arma::mat);
arma::mat eigvec_sym(arma::mat);
arma::mat inv_sympd(arma::mat);

#endif
