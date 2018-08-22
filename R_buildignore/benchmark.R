library(spectralGraphTopology)
library(ggplot2)

warmup_benchmark <- function(N_realizations, Nmax = 10) {
  T <- 50
  pb <- txtProgressBar(min = 3, max = Nmax, style = 3)
  N_set <- seq(3, Nmax, by = 1)
  rel_err_spec <- array(0, length(N_set))
  rel_err_naive <- array(0, length(N_set))
  rel_err_qp <- array(0, length(N_set))
  i <- 1
  # loop over the number of nodes
  for (N in N_set) {
    setTxtProgressBar(pb, N)
    rel_err <- c(naive = 0, qp = 0, spectral = 0)
    # loop over the number of realizations per node
    for (n in 1:N_realizations) {
      w <- runif(as.integer(.5 * N * (N - 1)))
      Lw <- L(w)
      Y <- MASS::mvrnorm(T, rep(0, N), MASS::ginv(Lw))
      covY <- cov(Y)
      res <- learnGraphTopology(Y, K = 1, beta = 50)
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
  x <- T / N_set
  df <- data.frame(x, rel_err_spec, rel_err_naive, rel_err_qp)
  ggplot(df, aes(x)) +
    geom_line(aes(y=rel_err_spec), colour="red") +
    geom_line(aes(y=rel_err_naive), colour="green") +
    geom_line(aes(y=rel_err_qp), colour="blue")
}

# naive estimator
naive <- function(Kmat) {
  return(MASS::ginv(Kmat))
}

# qp estimator
qp <- function(Kmat) {
  Sinv <- MASS::ginv(Kmat)
  R <- vecLmat(ncol(Sinv))
  qp <- quadprog::solve.QP(t(R) %*% R, t(R) %*% vec(Sinv), diag(ncol(R)))
  return (L(qp$solution))
}

# compute the relative error
relativeError <- function(Xtrue, Xest) {
    return (100 * norm(Xtrue - Xest, type = "F") / max(1., norm(Xtrue, type = "F")))
}

warmup_benchmark(N_realizations = 100, Nmax = 10)
warnings()
