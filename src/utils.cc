#include "utils.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;


// [[Rcpp::export]]
Eigen::MatrixXd blockDiagCpp(const std::vector<Eigen::MatrixXd>& matrices) {
  const unsigned int n = matrices.size();
  Eigen::VectorXd sizes(n);

  unsigned int N = 0;
  unsigned int cols, rows;
  for (unsigned int k = 0; k < n; ++k) {
    cols = matrices[k].cols();
    rows = matrices[k].rows();

    if(cols != rows) {
      std::stringstream err_msg;
      err_msg << "matrix is not square";
      throw std::invalid_argument(err_msg.str().c_str());
    } else {
      sizes(k) = cols;
      N += sizes(k);
    }
  }

  Eigen::MatrixXd blockDiag = Eigen::MatrixXd::Zero(N, N);
  unsigned int i = 0;
  for (unsigned int k = 0; k < n; ++k) {
    blockDiag.block(i, i, sizes(k), sizes(k)) = matrices[k];
    i += sizes(k);
  }
  return blockDiag;
}


// [[Rcpp::export]]
std::vector<double> metrics(const Eigen::MatrixXd& Wtrue, const Eigen::MatrixXd& West,
                            const double eps) {
  std::vector<double> metrics_list;
  double fscore, recall, specificity, accuracy;
  bool isthere_edge, isthere_est_edge;
  double tp = 0, fp = 0, fn = 0, tn = 0;
  const unsigned int n = Wtrue.cols();
  for (unsigned int i = 0; i < (n-1); ++i)
    for (unsigned int j = i+1; j < n; ++j) {
      isthere_edge = std::abs(Wtrue(i, j)) > eps;
      isthere_est_edge = std::abs(West(i, j)) > eps;
      if (isthere_edge && isthere_est_edge)
        tp += 1;
      else if (!isthere_edge && isthere_est_edge)
        fp += 1;
      else if (isthere_edge && !isthere_est_edge)
        fn += 1;
      else
        tn += 1;
    }
  fscore = 2 * tp / (2 * tp + fn + fp);
  recall = tp / (tp + fn);
  specificity = tn / (tn + fp);
  accuracy = (tp + tn) / (tp + tn + fp + fn);
  metrics_list.push_back(fscore);
  metrics_list.push_back(recall);
  metrics_list.push_back(specificity);
  metrics_list.push_back(accuracy);
  return metrics_list;
}


// [[Rcpp::export]]
Eigen::MatrixXd pairwise_matrix_rownorm2(const Eigen::MatrixXd& M) {
  const unsigned int n = M.rows();
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
  for (unsigned int i = 0; i < n-1; ++i)
    for (unsigned int j = i+1; j < n; ++j)
      V(i, j) = (M.row(i) - M.row(j)).squaredNorm();
  return V.selfadjointView<Upper>();
}


// [[Rcpp::export]]
Eigen::VectorXd upper_view_vec(const Eigen::MatrixXd& M) {
  const unsigned int p = M.cols();
  unsigned int t = 0;
  Eigen::VectorXd v(int(.5 * p * (p - 1)));
  for (unsigned int i = 0; i < p-1; ++i)
    for (unsigned int j = i+1; j < p; ++j) {
      v(t) = M(i, j);
      ++t;
    }
  return v;
}

