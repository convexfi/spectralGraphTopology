#include "operators.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd CppLOp(const Eigen::VectorXd& w, int n) {
    int j, k;
    Eigen::MatrixXd Lw = Eigen::MatrixXd::Zero(n, n);
    k = w.size();

    for (int i = n-2; i > -1; --i) {
        j = n - i - 1;
        Lw.row(i).tail(j) = -w.head(k).tail(j);
        k -= j;
    }

    Eigen::MatrixXd LwT = Lw.transpose();
    Lw += LwT;
    Eigen::MatrixXd LwColSum = Lw.colwise().sum().asDiagonal();
    Lw -= LwColSum;
    return Lw;
}

// [[Rcpp::export]]
Eigen::VectorXd CppLStarOp(const Eigen::MatrixXd& Y) {
    int n = Y.cols();
    int k = .5 * n * (n - 1);
    int j = 0;
    int l = 1;
    Eigen::VectorXd LStarY = Eigen::VectorXd::Zero(k);

    for (int i = 0; i < k; ++i) {
        LStarY(i) = Y(j, j) + Y(l, l) - (Y(l, j) + Y(j, l));
        if (l == (n - 1)) {
            l = (++j) + 1;
        } else {
            ++l;
        }
    }

    return LStarY;
}
