library(igraph)
library(corrplot)
library(spectralGraphTopology)
library(pals)
library(extrafont)
library(latex2exp)
set.seed(0)

n <- 49
n_realizations <- 10
p <- 30 * n
k <- 7
k_set <- c(1:14)
fs <- array(0, length(k_set))
re <- array(0, length(k_set))
bic <- array(0, length(k_set))
P <- diag(1, k)
# K-component graph
mgraph <- sample_sbm(n, pref.matrix = P, block.sizes = c(rep(n / k, k)))
# Erdo-Renyi as noise model
pe <- .45
a <- .25
erdos_renyi <- erdos.renyi.game(n, pe)
for (l in c(1:n_realizations)) {
  E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
  Ltrue <- as.matrix(laplacian_matrix(mgraph))
  E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
  Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
  Lnoisy <- Ltrue + Lerdo
  Y <- MASS::mvrnorm(p, mu = rep(0, n), Sigma = MASS::ginv(Lnoisy))
  S <- cov(Y)
  Sinv <- MASS::ginv(S)
  R <- spectralGraphTopology:::vecLmat(ncol(Sinv))
  qp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Sinv), diag(ncol(R)))
  w0 <- qp$solution
  for (j in c(1:length(k_set))) {
    cat("\nRunning simulation for K = ", k_set[j], "\n")
    graph <- learn_laplacian_matrix(S, k = k_set[j], w0 = w0, beta = 1e-1, alpha = 1e-2)
    re[j] <- re[j] + relativeError(Ltrue, graph$Lw)
    fs[j] <- fs[j] + Fscore(Ltrue, graph$Lw, 1e-2)
    bic[j] <- bic[j] + k_set[j] + graph$obj_fun[length(graph$obj_fun)]
  }
}
re <- re / n_realizations
fs <- fs / n_realizations
bic <- bic / n_realizations

saveRDS(re, "relerror.rds")
saveRDS(fs, "fscore.rds")
saveRDS(bic, "bic.rds")

