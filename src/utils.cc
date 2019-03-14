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


// TODO: put your sh*t together and write a class for this thing
// [[Rcpp::export]]
std::vector<double> metrics(const Eigen::MatrixXd& Wtrue, const Eigen::MatrixXd& West,
                            const double eps) {
  std::vector<double> metrics_list;
  double fscore, recall, specificity, accuracy;
  bool isthere_edge, isthere_est_edge;
  double tp = 0, fp = 0, fn = 0, tn = 0;
  const int n = Wtrue.cols();
  for (int i = 0; i < (n-1); ++i)
    for (int j = i+1; j < n; ++j) {
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
Eigen::MatrixXd pairwise_matrix_rownorm(const Eigen::MatrixXd& M) {
  const unsigned int n = M.rows();
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n-1; ++i)
    for (int j = i+1; j < n; ++j)
      V(i, j) = (M.row(i) - M.row(j)).squaredNorm();
  return V.selfadjointView<Upper>();
}
