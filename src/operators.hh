#ifndef OPERATORS_H
#define OPERATORS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;

MatrixXd L(const VectorXd&, int);
VectorXd Linv(const MatrixXd&);
VectorXd Lstar(const MatrixXd&);
VectorXd altLstar(const MatrixXd&);
MatrixXd vec(const MatrixXd&);
MatrixXd vecLmat(int);
#endif
