#ifndef UTILS_H
#define UTILS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Eigen;
using namespace std;

MatrixXd blockDiagCpp(const std::vector<MatrixXd>&);
std::vector<double> metrics(const MatrixXd&, const MatrixXd&, const double);
MatrixXd pairwise_matrix_rownorm(const MatrixXd&);
#endif
