#ifndef GAUSSIAN_PROCESS_H
#define GAUSSIAN_PROCESS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Eigen;

MatrixXd wiener_kernel(const unsigned int);
MatrixXd gaussian_kernel(const double, const double, const unsigned int);
MatrixXd rational_quadratic_kernel(const double, const double, const double,
                                   const unsigned int);
MatrixXd periodic_kernel(const double, const double, const double,
                         const unsigned int);
#endif
