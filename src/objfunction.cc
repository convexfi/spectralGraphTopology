#include "objfunction.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;


//' @export
// [[Rcpp::export]]
double objfunc(const Eigen::MatrixXd& Lw, const Eigen::MatrixXd& U,
               const Eigen::VectorXd& lambda, const Eigen::MatrixXd& Kmat,
               const double beta) {
  return loglikelihood(Lw, lambda, Kmat) + logprior(beta, Lw, lambda, U);
}

//' @export
// [[Rcpp::export]]
double loglikelihood(const Eigen::MatrixXd& Lw, const Eigen::VectorXd& lambda,
                     const Eigen::MatrixXd& Kmat) {
  return (Kmat * Lw).trace() - (lambda.array().log()).sum();
}

//' @export
// [[Rcpp::export]]
double logprior(const double beta, const Eigen::MatrixXd& Lw,
                const Eigen::VectorXd& lambda, const Eigen::MatrixXd& U) {
  Eigen::MatrixXd Uscaled = U * (lambda.array().sqrt()).matrix().asDiagonal();
  double frobnorm = (Lw -  Uscaled * Uscaled.transpose()).norm();
  return .5 * beta * frobnorm * frobnorm;
}
