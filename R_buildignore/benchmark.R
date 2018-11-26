# Run this file as Rscript benchmark.R

library(spectralGraphTopology)
library(viridis)
library(igraph)

# naive estimator
naive <- function(S) {
  return(MASS::ginv(S))
}

# qp estimator
qp <- function(S) {
  Sinv <- MASS::ginv(S)
  R <- vecLmat(ncol(Sinv))
  qp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Sinv), diag(ncol(R)))
  return (L(qp$solution))
}

# This function is a little benchmark to warm up
# against simpler estimators such as the naive one
# (generalized inverse of the covariance matrix)
# and one designed by Prof Dan which is the solution
# of a QP. Latter we will test against glasso and
# then against the (to be implemented) state-of-art
# algorithms.
warmup_benchmark <- function(N_realizations, N, ratios, beta_set) {
  P <- diag(.2 - .01, 4) + .01
  # fix the number of samples per node
  T_set <- as.integer(N * ratios)
  rel_err_spec <- matrix(0, ncol = length(T_set), nrow = length(beta_set))
  rel_err_naive <- array(0, length(T_set))
  rel_err_qp <- array(0, length(T_set))
  # loop over the number of nodes
  pb <- txtProgressBar(min = 1, max = N_realizations, style = 3)
  for (i in c(1:length(beta_set))) {
    for (j in c(1:length(T_set))) {
      cat("\nRunning simulation for", T_set[j], "samples per node, T/N = ", ratios[j],
          "and beta = ", beta_set[i], "\n")
      rel_err <- c(spectral = 0)
      # loop over the number of realizations per node
      for (n in 1:N_realizations) {
        setTxtProgressBar(pb, n)
        mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N/4, 4)))
        E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
        Lw <- as.matrix(laplacian_matrix(mgraph))
        Y <- MASS::mvrnorm(T_set[j], rep(0, N), MASS::ginv(Lw))
        #w <- runif(as.integer(.5 * N * (N - 1)))
        #Lw <- L(w)
        #Y <- MASS::mvrnorm(T_set[j], rep(0, N), MASS::ginv(Lw))
        #erdos_renyi <- erdos.renyi.game(N, .1)
        #E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 1e-1, max = 3)
        #Lw <- as.matrix(laplacian_matrix(erdos_renyi))
        #Y <- MASS::mvrnorm(n = T_set[j], mu = rep(0, N), Sigma = MASS::ginv(Lw))
        covY <- cov(Y)
        res <- learnGraphTopology(covY, K = 4, w0 = 'naive',
                                  beta = beta_set[i], maxiter = 50000)
        Lw_est <- res$Lw
        print(res$convergence)
        rel_err["spectral"] <- rel_err["spectral"] + relativeError(Lw, Lw_est)
        #rel_err["naive"] <- rel_err["naive"] + relativeError(Lw, naive(covY))
        #rel_err["qp"] <- rel_err["qp"] + relativeError(Lw, qp(covY))
      }
      # save the average of the relative errors
      rel_err_spec[i, j] <- rel_err["spectral"] / N_realizations
      #rel_err_naive[i] <- rel_err["naive"] / N_realizations
      #rel_err_qp[i] <- rel_err["qp"] / N_realizations
    }
  }
  colors = viridis(length(beta_set))
  plot(ratios, rel_err_spec[1, ], type = "b", pch=19, cex=.6,
       ylim=c(min(rel_err_spec) - 1e-2, max(rel_err_spec) + .1e-2),
       xlab = "T/N", ylab = "Relative Error", col = colors[1])
  for (i in c(2:length(beta_set))) {
    lines(ratios, rel_err_spec[i, ], type = "b", pch=19, cex=.6, col = colors[i])
  }
  legend("topright", legend=c("beta = .5",
                              "beta = 1",
                              "beta = 2.5",
                              "beta = 10",
                              "beta = 100"),
         col=colors, lty=c(1, 1, 1), cex=0.8)
}

# usage
ratios <- c(5, 20, 40, 60, 80, 100)
beta_set <- c(.5, 1., 2.5, 10, 100)
warmup_benchmark(N_realizations = 10, N = 64, ratios = ratios, beta_set = beta_set)
warnings()
