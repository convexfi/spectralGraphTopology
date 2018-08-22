library(spectralGraphTopology)
library(ggplot2)
set.seed(123)


# This function is a little benchmark to warm up
# against simpler estimators such as the naive one
# (generalized inverse of the covariance matrix)
# and one designed by Prof Dan which is the solution
# of a QP. Latter we will test against glasso and
# then against the (to be implemented) state-of-art
# algorithms.
warmup_benchmark <- function(N_realizations, T, ratios) {
  # fix the number of samples per node
  N_nodes_set <- sort(as.integer(T / ratios))
  rel_err_spec <- array(0, length(N_nodes_set))
  rel_err_naive <- array(0, length(N_nodes_set))
  rel_err_qp <- array(0, length(N_nodes_set))
  i <- 1
  # loop over the number of nodes
  pb <- txtProgressBar(min = 1, max = N_realizations, style = 3)
  for (N in N_nodes_set) {
    cat("\nRunning simultation for", N, "nodes T/N = ", ratios[length(N_nodes_set) - i + 1], "\n")
    rel_err <- c(naive = 0, qp = 0, spectral = 0)
    # loop over the number of realizations per node
    for (n in 1:N_realizations) {
      setTxtProgressBar(pb, n)
      w <- runif(as.integer(.5 * N * (N - 1)))
      Lw <- L(w)
      Y <- MASS::mvrnorm(T, rep(0, N), MASS::ginv(Lw))
      covY <- cov(Y)
      res <- learnGraphTopology(Y, K = 1, beta = 1, rho = .05,
                                maxiter_beta = 50)
      Lw_est <- res$Lw
      rel_err["spectral"] <- rel_err["spectral"] + relativeError(Lw, Lw_est)
      rel_err["naive"] <- rel_err["naive"] + relativeError(Lw, naive(covY))
      rel_err["qp"] <- rel_err["qp"] + relativeError(Lw, qp(covY))
    }
    # save the average of the relative errors
    rel_err_spec[i] <- rel_err["spectral"] / N_realizations
    rel_err_naive[i] <- rel_err["naive"] / N_realizations
    rel_err_qp[i] <- rel_err["qp"] / N_realizations
    i <- i + 1
  }
  df <- data.frame(ratios, rel_err_spec, rel_err_naive, rel_err_qp)
  ggplot(df, aes(ratios)) +
    geom_line(aes(y=rel_err_spec), colour="red") + geom_point(aes(y=rel_err_spec), colour="red") +
    geom_line(aes(y=rel_err_naive), colour="green") + geom_point(aes(y=rel_err_naive), colour="green") +
    geom_line(aes(y=rel_err_qp), colour="blue") + geom_point(aes(y=rel_err_qp), colour="blue") +
    labs(y="Relative Error (%)", x="T/N")
}

# usage
ratios <- c(.5, 1, 2, 5, 10, 20, 50)
warmup_benchmark(N_realizations = 100, T = 200, ratios = ratios)
warnings()
