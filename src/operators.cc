#include "operators.h"
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
  const int n = .5 * (1 + sqrt(1 + 8 * k));
  Eigen::MatrixXd Lw = Eigen::MatrixXd::Zero(n, n);

  for (int i = n-2; i > -1; --i) {
    j = n - i - 1;
    Lw.row(i).tail(j) = -w.head(k).tail(j);
    k -= j;
  }

  Eigen::MatrixXd LwColSum = (Lw + Lw.transpose()).colwise().sum().asDiagonal();
  Lw -= LwColSum;
  return Lw.selfadjointView<Upper>();
}

//' Computes the Adjacency linear operator which maps a vector of weights into
//' a valid Adjacency matrix.
//'
//' @param w weight vector of the graph
//' @return Aw the Adjacency matrix
//'
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd A(const Eigen::VectorXd& w) {
  int j;
  int k = w.size();
  const int N = .5 * (1 + sqrt(1 + 8 * k));
  Eigen::MatrixXd Aw = Eigen::MatrixXd::Zero(N, N);

  for (int i = N-2; i > -1; --i) {
    j = N - i - 1;
    Aw.row(i).tail(j) = w.head(k).tail(j);
    k -= j;
  }
  return Aw.selfadjointView<Upper>();
}


//' Computes the matrix form of the composition of the operators Lstar and
//' L, i.e., Lstar o L.
//'
//' @param n number of columns/rows
//' @return M the composition of Lstar and L
//'
// [[Rcpp::export]]
Eigen::MatrixXd Mmat(const int n) {
  Eigen::VectorXd e = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd M(n, n);
  e(0) = 1;
  M.col(0) = Lstar(L(e));
  for (int j = 1; j < n; ++j) {
    e(j - 1) = 0;
    e(j) = 1;
    M.col(j) = Lstar(L(e));
  }
  return M;
}


//' Computes the matrix form of the composition of the operators Astar and
//' A, i.e., Astar o A.
//'
//' @param n number of columns/rows
//' @return M the composition of Astar and A
//'
// [[Rcpp::export]]
Eigen::MatrixXd Pmat(const int n) {
  Eigen::VectorXd e = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd M(n, n);
  e(0) = 1;
  M.col(0) = Astar(A(e));
  for (int j = 1; j < n; ++j) {
    e(j - 1) = 0;
    e(j) = 1;
    M.col(j) = Astar(A(e));
  }
  return M;
}


//' Computes the matrix that represents the composition of
//' the vec and the L operators.
//'
//' @param n the dimension of L
//' @return R matrix such that vec(L(w)) = Rw
//'
// [[Rcpp::export]]
Eigen::MatrixXd vecLmat(const int n) {
  const int ncols = .5 * n * (n - 1);
  const int nrows = n * n;

  Eigen::VectorXd e = Eigen::VectorXd::Zero(ncols);
  Eigen::MatrixXd R(nrows, ncols);
  e(0) = 1;
  R.col(0) = vec(L(e));
  for (int j = 1; j < ncols; ++j) {
    e(j - 1) = 0;
    e(j) = 1;
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
// [[Rcpp::export]]
Eigen::VectorXd Lstar(const Eigen::MatrixXd& M) {
  int N = M.cols();
  int k = .5 * N * (N - 1);
  int j = 0;
  int l = 1;
  Eigen::VectorXd w(k);

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


//' Computes the Astar operator.
//'
//' @param M matrix
//' @return w vector
//'
// [[Rcpp::export]]
Eigen::VectorXd Astar(const Eigen::MatrixXd& M) {
  int N = M.cols();
  int k = .5 * N * (N - 1);
  int j = 0;
  int l = 1;
  Eigen::VectorXd w(k);

  for (int i = 0; i < k; ++i) {
    w(i) = M(l, j) + M(j, l);
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
// [[Rcpp::export]]
Eigen::VectorXd Linv(const Eigen::MatrixXd& M) {
  int N = M.cols();
  int k = .5 * N * (N - 1);
  Eigen::VectorXd w(k);
  int l = 0;

  for (int i = 0; i < N-1; ++i) {
    for (int j = i+1; j < N; ++j) {
      w(l) = -M(i, j);
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
// [[Rcpp::export]]
Eigen::VectorXd Ainv(const Eigen::MatrixXd& M) {
  int N = M.cols();
  int k = .5 * N * (N - 1);
  Eigen::VectorXd w(k);
  int l = 0;

  for (int i = 0; i < N-1; ++i) {
    for (int j = i+1; j < N; ++j) {
      w(l) = M(i, j);
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
  Eigen::VectorXd e = Eigen::VectorXd::Zero(k);

  e(0) = 1;
  w(0) = (MT * L(e)).trace();
  for (int i = 1; i < k; ++i) {
    e(i - 1) = 0;
    e(i) = 1;
    w(i) = (MT * L(e)).trace();
  }
  return w;
}

//' Alternative implementation for the Astar operator.
//' This is only used for unit testing. Use Astar for
//' a better performance.
//'
//' @param M matrix
//' @return w vector
//'
// [[Rcpp::export]]
Eigen::VectorXd altAstar(const Eigen::MatrixXd& M) {
  int N = M.cols();
  int k = .5 * N * (N - 1);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(k);
  Eigen::MatrixXd MT = M.transpose();
  Eigen::VectorXd e = Eigen::VectorXd::Zero(k);

  e(0) = 1;
  w(0) = (MT * A(e)).trace();
  for (int i = 1; i < k; ++i) {
    e(i - 1) = 0;
    e(i) = 1;
    w(i) = (MT * A(e)).trace();
  }
  return w;
}

