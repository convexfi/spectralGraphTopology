#include "utils.hh"
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
Eigen::MatrixXd blockDiagonal(const std::vector<Eigen::MatrixXd>& matrices) {
    int n = matrices.size();
    Eigen::VectorXd sizes(n);

    int N = 0;
    for (int k = 0; k < n; ++k) {
        sizes(k) = matrices[k].cols();
        N += sizes(k);
    }

    Eigen::MatrixXd blockDiag = Eigen::MatrixXd::Zero(N, N);
    int i = 0;
    for (int k = 0; k < n; ++k) {
        blockDiag.block(i, i, sizes(k), sizes(k)) = matrices[k];
        i += sizes(k);
    }
    return blockDiag;
}
