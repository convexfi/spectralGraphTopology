#ifndef ADMM_H
#define ADMM_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Eigen;
using namespace std;

Eigen::MatrixXd Theta_update(const Eigen::MatrixXd&, const Eigen::VectorXd&,
                             const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                             const double);
#endif
