#include "operators.hh"
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;


//' Computes the Laplacian linear operator which maps a vector of weights into
//' a valid Laplacian matrix.
//'
//' @param w weight vector of the graph
//' @return Lw the Laplacian matrix
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd L(const Eigen::VectorXd& w) {
    int j;
    int k = w.size();
    int N = .5 * (1 + sqrt(1 + 8 * k));
    Eigen::MatrixXd Lw = Eigen::MatrixXd::Zero(N, N);

    for (int i = N-2; i > -1; --i) {
        j = N - i - 1;
        Lw.row(i).tail(j) = -w.head(k).tail(j);
        k -= j;
    }

    Eigen::MatrixXd LwT = Lw.transpose();
    Lw += LwT;
    Eigen::MatrixXd LwColSum = Lw.colwise().sum().asDiagonal();
    Lw -= LwColSum;
    return Lw;
}


//' Computes the matrix that represents the composition of
//' the vec and the L operators.
//'
//' @param n the dimension of L
//' @return R matrix such that vec(L(w)) = Rw
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd vecLmat(int n) {
    int ncols = .5 * n * (n - 1);
    int nrows = n * n;

    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(nrows, ncols);
    for (int j = 0; j < ncols; ++j) {
        Eigen::VectorXd e = Eigen::VectorXd::Zero(ncols);
        e(j) = 1.;
        R.col(j) = vec(L(e));
    }
    return R;
}


//' Computes the canonical vec operator, i.e., converts a given n x m matrix M,
//' into a  nm x 1 column vector vec(M) by stacking the columns of M on top of
//' one another.
//'
//' @param M input matrix
//' @return w vector such that w = vec(M)
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd vec(const Eigen::MatrixXd& M) {
    Eigen::MatrixXd w(M.size(), 1);
    int k = 0;
    for (int j = 0; j < M.cols(); ++j) {
        for (int i = 0; i < M.rows(); ++i) {
            w(k, 0) = M(i, j);
            ++k;
        }
    }
    return w;
}


//' Computes the Lstar operator.
//'
//' @param M matrix
//' @return w vector
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd Lstar(const Eigen::MatrixXd& M) {
    int N = M.cols();
    int k = .5 * N * (N - 1);
    int j = 0;
    int l = 1;
    Eigen::VectorXd w = Eigen::VectorXd::Zero(k);

    for (int i = 0; i < k; ++i) {
        w(i) = M(j, j) + M(l, l) - (M(l, j) + M(j, l));
        if (l == (N - 1)) {
            l = (++j) + 1;
        } else {
            ++l;
        }
    }
    return w;
}


//' Computes the inverse of the L operator.
//'
//' @param M Laplacian matrix
//' @return w the weight vector of the graph
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd Linv(const Eigen::MatrixXd& M) {
    int N = M.cols();
    int k = .5 * N * (N - 1);
    Eigen::VectorXd w = Eigen::VectorXd::Zero(k);
    int l = 0;

    for (int i = 0; i < N-1; ++i) {
        for (int j = i+1; j < N; ++j) {
            w(l) = -M(i, j);
            ++l;
        }
    }
    return w;
}


//' Alternative implementation for the Lstar operator.
//' This is only used for unit testing. Use Lstar for
//' a better performance.
//'
//' @param M matrix
//' @return w vector
//'
// [[Rcpp::export]]
Eigen::VectorXd altLstar(const Eigen::MatrixXd& M) {
    int N = M.cols();
    int k = .5 * N * (N - 1);
    Eigen::VectorXd w = Eigen::VectorXd::Zero(k);
    Eigen::MatrixXd MT = M.transpose();

    for (int i = 0; i < k; ++i) {
        Eigen::VectorXd e = Eigen::VectorXd::Zero(k);
        e(i) = 1.;
        w(i) = (MT * L(e)).trace();
    }
    return w;
}
