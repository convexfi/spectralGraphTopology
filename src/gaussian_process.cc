#define _USE_MATH_DEFINES
#include <cmath>
#include "gaussian_process.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;


//' @export
// [[Rcpp::export]]
Eigen::MatrixXd wiener_kernel(const unsigned int n) {
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = i; j < n; ++j)
      K(i, j) = std::min(i+1, j+1);
  return K.selfadjointView<Upper>();
}

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd gaussian_kernel(const double a, const double s,
                                const unsigned int n) {
  // a and s have to be positive
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
  int d;
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = i; j < n; ++j) {
      d = i - j;
      K(i, j) = a * std::exp(-s * d * d);
    }
  return K.selfadjointView<Upper>();
}


//' @export
// [[Rcpp::export]]
Eigen::MatrixXd rational_quadratic_kernel(const double a, const double b,
                                          const double s, const unsigned int n) {
  // a, b, and s have to be positive
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
  int d;
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = i; j < n; ++j) {
      d = i - j;
      K(i, j) = a * std::pow(1 + s * d * d / b, -b);
    }
  return K.selfadjointView<Upper>();
}


//' @export
// [[Rcpp::export]]
Eigen::MatrixXd periodic_kernel(const double a, const double p,
                                const double s, const unsigned int n) {
  // a, s have to be positive
  double sin_t;
  int d;
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n, n);
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = i; j < n; ++j) {
      d = i - j;
      sin_t = std::sin(M_PI * std::abs(d) / p);
      K(i, j) = a * std::exp(-s * sin_t * sin_t);
    }
  return K.selfadjointView<Upper>();
}
