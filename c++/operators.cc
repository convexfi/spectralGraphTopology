#include <operators.hh>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd LOp(VectorXd w, int n) {
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

// [[Rcpp::export]]
VectorXd LStarOp(MatrixXd Y) {
    int n = Y.cols();
    int k = .5 * n * (n - 1);
    VectorXd LStarY = VectorXd::Zero(k);
    MatrixXd Lw = MatrixXd::Zero(n, n);

    for (int i = 0; i < k; ++i) {
        VectorXd w = VectorXd::Zero(k);
        w(i) = 1.;
        Lw = LOp(w, n);
        LStarY(i) = (Y.transpose() * Lw).trace();
    }

    return LStarY;
}
