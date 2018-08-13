#ifndef OPERATORS_H
#define OPERATORS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Eigen;

MatrixXd CppLOp(const Eigen::VectorXd&, int);
VectorXd CppLStarOp(const Eigen::MatrixXd&);
VectorXd CppLStarOpImpl(const Eigen::MatrixXd&);
#endif
