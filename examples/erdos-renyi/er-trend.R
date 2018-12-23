library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(0)
N_realizations <- 10
ratio <- 100
N <- 32
T <- ratio * N
p <- .15
erdos_renyi <- erdos.renyi.game(N, p)
Ltrue <- as.matrix(laplacian_matrix(erdos_renyi))
# sample data from GP with covariance matrix set as
# the pseudo inverse of the true Laplacian
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
S <- cov(Y)
# run spectralGraphTopology
graph <- learnGraphTopology(S, w0 = "qp", beta = 10, maxiter = 100000, ftol = 1e-8, Lwtol = 1e-8)
print(graph$convergence)
print(relativeError(Ltrue, graph$Lw))
print(Fscore(Ltrue, graph$Lw, 1e-1))

gr = .5 * (1 + sqrt(5))
colors <- c("#706FD3", "#FF5252", "#33D9B2")
setEPS()
postscript("er_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:length(graph$loglike)), graph$loglike, type = "b", lty = 1, pch = 15, cex=.75, col = colors[1],
     xlab = "Iteration Number", ylab = "")
grid()
lines(c(1:length(graph$loglike)), graph$obj_fun, type = "b", xaxt = "n", lty = 2, pch=16, cex=.75, col = colors[2])
lines(c(1:length(graph$loglike)), graph$obj_fun - graph$loglike, type = "b", xaxt = "n", lty = 3, pch=17, cex=.75,
      col = colors[3])
legend("topright", legend = c("likelihood", "posterior", "prior"),
       col=colors, pch=c(15, 16, 17), lty=c(1, 2, 3), bty="n")
dev.off()
embed_fonts("er_trend.ps", outfile="er_trend.ps")
