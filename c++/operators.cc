#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;

// [[Rcpp::export]]
MatrixXd LOp(Map<VectorXd> w, int n) {
    int j, k;
    MatrixXd Lw = MatrixXd::Zero(n, n);
    k = w.size();
    for (int i = n-2; i > -1; --i) {
        j = n - i - 1;
        Lw.row(i).tail(j) = -w.head(k).tail(j);
        k -= j;
    }
    MatrixXd LwT = Lw.transpose();
    Lw += LwT;
    MatrixXd LwColSum = Lw.colwise().sum().asDiagonal();
    Lw -= LwColSum;
    return Lw;
}
