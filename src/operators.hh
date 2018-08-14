#ifndef OPERATORS_H
#define OPERATORS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;

MatrixXd L(const Eigen::VectorXd&, int);
VectorXd Lstar(const Eigen::MatrixXd&);
VectorXd altLstar(const Eigen::MatrixXd&);
#endif
