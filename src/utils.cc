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
