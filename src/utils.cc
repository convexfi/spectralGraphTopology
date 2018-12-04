#include "utils.hh"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;


// [[Rcpp::export]]
Eigen::MatrixXd blockDiagCpp(const std::vector<Eigen::MatrixXd>& matrices) {
  const int n = matrices.size();
  Eigen::VectorXd sizes(n);

  int N = 0;
  int cols, rows;
  for (int k = 0; k < n; ++k) {
    cols = matrices[k].cols();
    rows = matrices[k].rows();

    if(cols != rows) {
      std::stringstream err_msg;
      err_msg << "matrix " << std::to_string(k + 1) << " is not square";
      throw std::invalid_argument(err_msg.str().c_str());
    } else {
      sizes(k) = cols;
      N += sizes(k);
    }
  }

  Eigen::MatrixXd blockDiag = Eigen::MatrixXd::Zero(N, N);
  int i = 0;
  for (int k = 0; k < n; ++k) {
    blockDiag.block(i, i, sizes(k), sizes(k)) = matrices[k];
    i += sizes(k);
  }
  return blockDiag;
}

// [[Rcpp::export]]
double Fscore(const Eigen::MatrixXd& Wtrue, const Eigen::MatrixXd& West,
              const double eps) {
  bool isthere_edge, isthere_est_edge;
  double tp = 0, fp = 0, fn = 0;
  const int n = Wtrue.cols();
  for (int i = 0; i < (n-1); ++i) {
    for (int j = i+1; j < n; ++j) {
      isthere_edge = std::abs(Wtrue(i, j)) > eps;
      isthere_est_edge = std::abs(West(i, j)) > eps;
      if (isthere_edge && isthere_est_edge)
        tp += 1;
      else if (!isthere_edge && isthere_est_edge)
        fp += 1;
      else if (isthere_edge && !isthere_est_edge)
        fn += 1;
    }
  }

  return 2 * tp / (2 * tp + fn + fp);
}
