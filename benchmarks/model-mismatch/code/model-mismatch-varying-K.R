library(igraph)
library(corrplot)
library(spectralGraphTopology)
library(pals)
library(extrafont)
library(latex2exp)
set.seed(0)

N <- 49
Nrealizations <- 10
T <- 30 * N
K <- 7
K_set <- c(1:7)
fs <- array(0, length(K_set))
re <- array(0, length(K_set))
P <- diag(1, K)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
# Erdo-Renyi as noise model
p <- .45
a <- .25
erdos_renyi <- erdos.renyi.game(N, p)
for (n in c(1:Nrealizations)) {
  E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
  Ltrue <- as.matrix(laplacian_matrix(mgraph))
  E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
  Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
  Lnoisy <- Ltrue + Lerdo
  Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Lnoisy))
  S <- cov(Y)
  Sinv <- MASS::ginv(S)
  R <- vecLmat(ncol(Sinv))
  qp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Sinv), diag(ncol(R)))
  w0 <- qp$solution
  for (j in c(1:length(K_set))) {
    cat("\nRunning simulation for K = ", K_set[j], "\n")
    graph <- learnGraphTopology(S, K = K_set[j], w0 = w0, beta = 1e-1,
                                alpha = 1e-2, ftol = 1e-4, Lwtol = 1e-4)
    re[j] <- re[j] + relativeError(Ltrue, graph$Lw)
    fs[j] <- fs[j] + Fscore(Ltrue, graph$Lw, 1e-2)
  }
}
re <- re / Nrealizations
fs <- fs / Nrealizations
print(re)

gr = .5 * (1 + sqrt(5))
colors <- c("#706FD3")
setEPS()
postscript("../latex/figures/fscore.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(K_set, fs, type = "b", pch=15, cex=.75, ylim=c(min(fs), max(fs)),
     xlab = TeX("$\\mathit{K}$"), ylab = "Average F-score", col = colors[1])
grid()
legend("topleft", legend = c("SGL"), col=c(colors[1]), pch=c(15), lty=c(1), bty="n")
dev.off()
embed_fonts("../latex/figures/fscore.ps", outfile="../latex/figures/fscore.ps")
