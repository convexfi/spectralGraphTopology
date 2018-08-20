#ifndef UTILS_H
#define UTILS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace std;

MatrixXd blockDiag(const std::vector<MatrixXd>&);
#endif
