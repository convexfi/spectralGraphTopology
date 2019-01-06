#ifndef OPERATORS_H
#define OPERATORS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;

// Laplacian matrix operators
MatrixXd L(const VectorXd&);
VectorXd Linv(const MatrixXd&);
VectorXd Lstar(const MatrixXd&);
MatrixXd Mmat(const int);
VectorXd altLstar(const MatrixXd&);
MatrixXd vecLmat(const int);
// Adjacency matrix operators
MatrixXd A(const VectorXd&);
VectorXd Ainv(const MatrixXd&);
VectorXd Astar(const MatrixXd&);
MatrixXd Pmat(const int);
VectorXd altAstar(const MatrixXd&);
// Utility operators
MatrixXd vec(const MatrixXd&);
#endif
