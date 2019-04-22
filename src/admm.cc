#include "admm.hh"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd Theta_update(const Eigen::MatrixXd& U,
                             const Eigen::VectorXd& lambda,
                             const Eigen::MatrixXd& K,
                             const Eigen::MatrixXd& M,
                             const double beta) {
  const int p = U.rows();
  const int v = U.cols();
  Eigen::MatrixXd lmdUT(v, p);
  Eigen::MatrixXd Ulmd(p, v);
  Eigen::MatrixXd X(p, p);
  Eigen::MatrixXd Theta = Eigen::MatrixXd::Ones(p, p);

  for (int i = 0; i < v; ++i)
    lmdUT.row(i) = U.col(i) * std::sqrt(lambda(i));
  Ulmd = lmdUT.transpose();
  X = Ulmd * lmdUT - (M + K) / beta;

  for (int i = 0; i < (p-1); ++i)
    for (int j = i+1; j < p; ++j)
      Theta(i, j) = std::min(std::max(X(i, j), -1.0), 0.0);

  return Theta.selfadjointView<Upper>();
}
