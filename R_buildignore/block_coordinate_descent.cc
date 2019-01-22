#include "block_coordinate_descent.hh"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

Eigen::VectorXd w_update_cpp(const Eigen::VectorXd& w, const Eigen::MatrixXd& Lw,
                             const Eigen::MatrixXd& U, const double beta,
                             const Eigen::VectorXd& lambda,
                             const Eigen::MatrixXd& Kmat) {
  int k = w.length();
  int Ucols = U.cols();
  Eigen::VectorXd w_update(k);
  Eigen::MatrixXd lmdUT(Ucols, U.rows());

  for (int i = 0; i < Ucols; ++i)
    lmdUT.row(i) = U.col(i) * std::sqrt(lambda(i));

  w_update = w - (.5 / Ucols) * Lstar(Lw - lmdUT.transpose() * lmdUT + Kmat / beta);
  for (int i = 0; i < k; ++i)
    if (w_update(i) < 0)
      w_update(i) = 0;
  return w_udpate;
}


Eigen::VectorXd lambda_update_cpp(const Eigen::MatrixXd& U,
                                  const Eigen::MatrixXd& Lw,
                                  const Eigen::MatrixXd& Lw,
                                  const double beta, const unsigned int n, const unsigned int k) {
  const unsigned int q = n - k;
  Eigen::ArrayXd d(q), lmd(q);
  d = (U.transpose() * Lw * U).diagonal().array();
  lmd = .5 * (d - std::sqrt(d * d + 4 / beta));

  for (unsigned int i = 0; i < q; ++i) {
    if (lmd(i) < lb)
      lmd(i) = lb;
    if (lmd(i) > ub)
  }
}
