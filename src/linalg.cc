#include <RcppArmadillo.h>
#include "linalg.hh"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::vec eigenvalues(arma::mat M) {
  return arma::eig_sym(M);
}

// [[Rcpp::export]]
arma::mat eigenvectors(arma::mat M) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, M);
  return eigvec;
}

// [[Rcpp::export]]
arma::mat inv_sympd(arma::mat M) {
  return arma::inv_sympd(M);
}
