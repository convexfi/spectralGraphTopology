#ifndef OPERATORS_H
#define OPERATORS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;

MatrixXd L(const VectorXd&);
VectorXd Linv(const MatrixXd&);
VectorXd Lstar(const MatrixXd&);
MatrixXd Mmat(const int);
VectorXd altLstar(const MatrixXd&);
MatrixXd vec(const MatrixXd&);
MatrixXd vecLmat(const int);
MatrixXd A(const VectorXd&);
#endif
