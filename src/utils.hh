#ifndef UTILS_H
#define UTILS_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace std;

MatrixXd blockDiagCpp(const std::vector<MatrixXd>&);
double Fscore(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const double);
#endif
