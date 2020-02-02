#' @title Cluster a k-component graph from data using the Constrained Laplacian Rank algorithm
#'
#' Cluster a k-component graph on the basis of an observed data matrix.
#' Check out https://mirca.github.io/spectralGraphTopology for code examples.
#'
#' @param Y a pxn data matrix, where p is the number of nodes and n is the number of
#'        features (or data points per node)
#' @param k the number of components of the graph
#' @param m the maximum number of possible connections for a given node used
#'        to build an affinity matrix
#' @param lmd L2-norm regularization hyperparameter
#' @param eigtol value below which eigenvalues are considered to be zero
#' @param edgetol value below which edge weights are considered to be zero
#' @param maxiter the maximum number of iterations
#' @return A list containing the following elements:
#' \item{\code{Laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{Adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{eigvals}}{the eigenvalues of the Laplacian Matrix}
#' \item{\code{lmd_seq}}{sequence of lmd values at every iteration}
#' \item{\code{elapsed_time}}{elapsed time at every iteration}
#' @author Ze Vinicius and Daniel Palomar
#' @references Nie, Feiping and Wang, Xiaoqian and Jordan, Michael I. and Huang, Heng.
#'             The Constrained Laplacian Rank Algorithm for Graph-based Clustering, 2016,
#'             AAAI'16. http://dl.acm.org/citation.cfm?id=3016100.3016174
#' @examples
#' library(clusterSim)
#' library(spectralGraphTopology)
#' library(igraph)
#' set.seed(1)
#' # number of nodes per cluster
#' N <- 30
#' # generate datapoints
#' twomoon <- shapes.two.moon(N)
#' # estimate underlying graph
#' graph <- cluster_k_component_graph(twomoon$data, k = 2)
#' # build network
#' net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
#' # colorify nodes and edges
#' colors <- c("#706FD3", "#FF5252", "#33D9B2")
#' V(net)$cluster <- twomoon$clusters
#' E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
#'                       function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
#'                                         colors[V(net)$cluster[x[1]]], '#000000'))
#' V(net)$color <- c(colors[1], colors[2])[twomoon$clusters]
#' # plot network
#' plot(net, layout = twomoon$data, vertex.label = NA, vertex.size = 3)
#' @export
cluster_k_component_graph <- function(Y, k = 1, m = 5, lmd = 1, eigtol = 1e-9,
                                      edgetol = 1e-6, maxiter = 1000) {
  time_seq <- c(0)
  start_time <- proc.time()[3]
  A <- build_initial_graph(Y, m)
  n <- ncol(A)
  S <- matrix(1/n, n, n)
  DS <- diag(.5 * colSums(S + t(S)))
  LS <-  DS - .5 * (S + t(S))
  DA <- diag(.5 * colSums(A + t(A)))
  LA <- DA - .5 * (A + t(A))
  if (k == 1)
    F <- matrix(eigvec_sym(LA)[, 1:k])
  else
    F <- eigvec_sym(LA)[, 1:k]
  # bounds for variables in the QP solver
  bvec <- c(1, rep(0, n))
  Amat <- cbind(rep(1, n), diag(n))
  lmd_seq <- c(lmd)
  pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  lambda: :lmd  null_eigvals: :null_eigvals",
                                   total = maxiter, clear = FALSE, width = 100)
  for (ii in c(1:maxiter)) {
    V <- pairwise_matrix_rownorm2(F)
    for (i in c(1:n)) {
      p <- A[i, ] - .5 * lmd * V[i, ]
      qp <- quadprog::solve.QP(Dmat = diag(n), dvec = p, Amat = Amat, bvec = bvec, meq = 1)
      S[i, ] <- qp$solution
    }
    DS <- diag(.5 * colSums(S + t(S)))
    LS <- DS - .5 * (S + t(S))
    F <- eigvec_sym(LS)[, 1:k]
    eig_vals <- eigval_sym(LS)
    n_zero_eigenvalues <- sum(abs(eig_vals) < eigtol)
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    pb$tick(token = list(lmd = lmd, null_eigvals = n_zero_eigenvalues))
    if (k < n_zero_eigenvalues)
      lmd <- .5 * lmd
    else if (k > n_zero_eigenvalues)
      lmd <- 2 * lmd
    else
      break
    lmd_seq <- c(lmd_seq, lmd)
  }
  LS[abs(LS) < edgetol] <- 0
  AS <- diag(diag(LS)) - LS
  return(list(Laplacian = LS, Adjacency = AS, eigenvalues = eig_vals,
              lmd_seq = lmd_seq, elapsed_time = time_seq))
}

build_initial_graph <- function(Y, m) {
  n <- nrow(Y)
  A <- matrix(0, n, n)
  E <- pairwise_matrix_rownorm2(Y)
  for (i in c(1:n)) {
    sorted_index <- order(E[i, ])
    j_sweep <- sorted_index[2:(m+1)]
    den <- m * E[i, sorted_index[m+2]] - sum(E[i, j_sweep])
    ei <- E[i, sorted_index[m+2]]
    for (j in j_sweep) {
      A[i, j] <- (ei - E[i, j]) / den
    }
  }
  return(A)
}
