#include "linalg.hh"

using namespace OsqpEigen;
using namespace arma;
using namespace Eigen;

Eigen::MatrixXd constr_laplacian_rank(const arma::mat& A, const unsigned int k,
                                      const double lmd,
                                      const double eig_tol = 1e-4,
                                      const unsigned int maxiter = 1000,
                                      const unsigned int type = 1) {
  const unsigned int n = A.n_cols();
  arma::mat D = arma::zeros<arma::mat>(n, n);
  Eigen::VectorXd p = Eigen::VectorXd(n);
  Eigen::MatrixXd V = Eigen::MatrixXd(n, n);
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(n, n);
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
  D.diagonal() = .5 * (arma::sum(A) + arma::sum(A, 1));
  arma::mat F = eigenvectors(D - .5 * (A.t() + A)).head_cols(k);

  // bounds for variables in the QP solver
  Eigen::VectorXd l = Eigen::VectorXd::Zero(n+1);
  Eigen::VectorXd u = Eigen::VectorXd::Infinity(n+1);
  Eigen::MatrixXd Amat = Eigen::MatrixXd(n+1, n);
  Amat << Eigen::MatrixXd::Identity(n, n), Eigen::VectorXd::One(n);
  l(n) = 1;
  u(n) = 1;

  // instantiate solver
  OsqpEigen::Solver solver;

  for (unsigned int ii = 0; ii < maxiter; ++ii) {
    // update S
    if (type == 1) {
      V = pairwise_matrix_row_norm(F);
      for (unsigned int i = 0; i < n; ++i) {
        U.diagonal() = 1 / (S.row(i) - A.row(i)).abs();
        p = .5 * lmd * V.row(i) - U * A.row(i);
        S.row(i) = solve_osqp(U, p, Amat, l, u);
      }
    } else if (type == 2) {
      for (unsigned int i = 0; i < n; ++i) {
        S.row(i) =
      }
    }
    arma::mat S_arma = arma::mat(S.data(), n, n, false, false);
    D.diag() = .5 * (arma::sum(S_arma) + arma::sum(S_arma, 1));
    // update F
    F = eigenvectors(D - .5 * (S_arma.t() + S_arma)).head_cols(k);
    if (arma::sum(eigenvalues(S_arma).head(k)) < k * eig_tol)
      break;
  }
  return S;
}

Eigen::MatrixXd pairwise_matrix_row_norm(const Eigen::MatrixXd& M) {
  const unsigned int n = M.rows();
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);

  for (unsigned int i = 0; i < (n-1); ++i)
    for (unsigned int j = (i+1); j < n; ++j)
      V(i, j) = (M.row(i) - M.row(j)).norm().square();

  return V.selfadjointView<Upper>();
}

Eigen::VectorXd solve_osqp(const Eigen::MatrixXd& P, const Eigen::VectorXd& q,
                           const Eigen::MatrixXd& A, const Eigen::VectorXd& l, const Eigen::VectorXd& u) {

}
