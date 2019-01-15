library(igraph)
library(viridis)
library(corrplot)
library(spectralGraphTopology)
library(latex2exp)
library(extrafont)
set.seed(0)

N_realizations <- 10
N <- 16
T <- 100 * N
K <- 4
P <- diag(1, K)
# probability of node connection
p <- .25
kappa_seq <- c(.1, .2, .3, .4, .5, .6)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
rel_err <- array(0, length(kappa_seq))
fscore <- array(0, length(kappa_seq))
for (k in c(1:length(kappa_seq))) {
  a <- kappa_seq[k]
  cat("\nRunning simulation for kappa =", a, "\n")
  for (n in 1:N_realizations) {
    E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Wtrue <- diag(diag(Ltrue)) - Ltrue
    # Erdo-Renyi as noise model
    erdos_renyi <- erdos.renyi.game(N, p)
    E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
    Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
    # Noisy Laplacian
    Lnoisy <- Ltrue + Lerdo
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Lnoisy))
    S <- cov(Y)

    Sinv <- MASS::ginv(S)
    R <- vecLmat(ncol(Sinv))
    qp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Sinv), diag(ncol(R)))
    w0 <- qp$solution

    graph <- learnLaplacianGraphTopology(S, w0 = w0, K = K, beta = 100 * N,
                                alpha = 1e-2, maxiter = 100000)
    print(graph$convergence)
    fs <- Fscore(Ltrue, graph$Lw, 1e-3)
    re = relativeError(Ltrue, graph$Lw)
    rel_err[k] <- rel_err[k] + re
    fscore[k] <- fscore[k] + fs
  }
  rel_err[k] <- rel_err[k] / N_realizations
  fscore[k] <- fscore[k] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err[k], "\n")
  cat("Avg Fscore: ")
  cat(fscore[k], "\n")
}
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("..latex/figures/relative_error_kappa.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:length(kappa_seq)), rel_err, type = "b", pch=19, cex=.75, ylim=c(min(rel_err), max(rel_err)),
     xlab = TeX("$\\kappa$"), ylab = "Average Relative Error", col = "#706FD3", xaxt = "n")
grid()
axis(side = 1, at = c(1:length(kappa_seq)), labels = kappa_seq)
dev.off()
embed_fonts("..latex/figures/relative_error_kappa.ps", outfile="..latex/figures/relative_error_kappa.ps")
setEPS()
postscript("..latex/figures/fscore_kappa.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:length(kappa_seq)), fscore, type = "b", pch=19, cex=.75, ylim=c(min(fscore), 1.),
     xlab = TeX("$\\kappa$"), ylab = "Average F-score", col = "#706FD3", xaxt = "n")
grid()
axis(side = 1, at = c(1:length(kappa_seq)), labels = kappa_seq)
dev.off()
embed_fonts("..latex/figures/fscore_kappa.ps", outfile="..latex/figures/fscore_kappa.ps")
