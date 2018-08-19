#ifndef UTILS_H
#define UTILS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace std;

MatrixXd blockDiagonal(const std::vector<MatrixXd>& matrices);
#endif
