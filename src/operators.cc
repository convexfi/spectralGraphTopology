#include "operators.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd L(const Eigen::VectorXd& w) {
   /* Computes the L Operator.
    *
    * Args:
    *   w: weight vector of the graph
    *
    * Returns:
    *   Theta: the Laplacian matrix
    **/
    int j;
    int k = w.size();
    int N = .5 * (1 + sqrt(1 + 8 * k));
    Eigen::MatrixXd Theta = Eigen::MatrixXd::Zero(N, N);

    for (int i = N-2; i > -1; --i) {
        j = N - i - 1;
        Theta.row(i).tail(j) = -w.head(k).tail(j);
        k -= j;
    }

    Eigen::MatrixXd ThetaT = Theta.transpose();
    Theta += ThetaT;
    Eigen::MatrixXd ThetaColSum = Theta.colwise().sum().asDiagonal();
    Theta -= ThetaColSum;
    return Theta;
}

// [[Rcpp::export]]
Eigen::VectorXd Lstar(const Eigen::MatrixXd& Y) {
   /* Computes the Lstar operator.
    *
    * Args:
    *   Y: matrix
    *
    * Returns:
    *   w: vector
    */
    int N = Y.cols();
    int k = .5 * N * (N - 1);
    int j = 0;
    int l = 1;
    Eigen::VectorXd w = Eigen::VectorXd::Zero(k);

    for (int i = 0; i < k; ++i) {
        w(i) = Y(j, j) + Y(l, l) - (Y(l, j) + Y(j, l));
        if (l == (N - 1)) {
            l = (++j) + 1;
        } else {
            ++l;
        }
    }
    return w;
}

// [[Rcpp::export]]
Eigen::VectorXd Linv(const Eigen::MatrixXd& Theta) {
   /* Computes the inverse of the L Operator.
    *
    * Args:
    *   Theta: the Laplacian matrix
    *
    * Returns:
    *   w: weight vector of the graph
    **/
    int N = Theta.cols();
    int k = .5 * N * (N - 1);
    Eigen::VectorXd w = Eigen::VectorXd::Zero(k);
    int l = 0;

    for (int i = 0; i < N-1; ++i) {
        for (int j = i+1; j < N; ++j) {
            w(l) = -Theta(i, j);
            ++l;
        }
    }

    return w;
}

// [[Rcpp::export]]
Eigen::VectorXd altLstar(const Eigen::MatrixXd& Y) {
   /* Alternative implementation for the Lstar operator.
    * This is only used for unit testing. Use Lstar for
    * a better performance.
    *
    *  Args:
    *    Y: matrix
    *
    *  Returns:
    *    w: vector
    */
    int N = Y.cols();
    int k = .5 * N * (N - 1);
    Eigen::VectorXd w = Eigen::VectorXd::Zero(k);
    Eigen::MatrixXd YT = Y.transpose();

    for (int i = 0; i < k; ++i) {
        Eigen::VectorXd e = Eigen::VectorXd::Zero(k);
        e(i) = 1.;
        w(i) = (YT * L(e)).trace();
    }
    return w;
}
