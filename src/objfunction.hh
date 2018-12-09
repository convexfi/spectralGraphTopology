#ifndef OBJFUNCTION_H
#define OBJFUNCTION_H
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Eigen;
double objfunc(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
               const Eigen::VectorXd&, const Eigen::MatrixXd&, const double beta);
double loglikelihood(const Eigen::MatrixXd&, const Eigen::VectorXd&,
                     const Eigen::MatrixXd&);
double logprior(const double, const Eigen::MatrixXd&, const Eigen::VectorXd&,
                const Eigen::MatrixXd&);
#endif
