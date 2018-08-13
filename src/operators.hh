#ifndef OPERATORS_H
#define OPERATORS_H
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;

MatrixXd CppLOp(VectorXd w, int n);
VectorXd CppLStarOp(MatrixXd Y);
#endif
