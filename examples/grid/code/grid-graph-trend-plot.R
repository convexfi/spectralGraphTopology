library(igraph)
library(spectralGraphTopology)
library(scales)
library(extrafont)

set.seed(0)
# number of nodes
N <- 64
# generate structure
grid <- make_lattice(length = sqrt(N), dim = 2)
# number of samples
T <- as.integer(100 * N)
# set edges weights randomly
E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
# get true Laplacian and true Adjacency matrix
Ltrue <- as.matrix(laplacian_matrix(grid))
Wtrue <- diag(diag(Ltrue)) - Ltrue
# sample data
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
# estimate underlying graph
graph <- learnGraphTopology(cov(Y), w0 = "naive", beta = 10, alpha = 5e-3, Lwtol = 1e-4)
# print convergence, relative error and fscore
cat("\nconvergence status: ", graph$convergence, "\n")
cat("relative error: ", relativeError(Ltrue, graph$Lw), "\n")
cat("fscore: ", Fscore(Ltrue, graph$Lw, 1e-1), "\n")
# plot convergence trend
gr = .5 * (1 + sqrt(5))
colors <- c("#706FD3", "#FF5252", "#33D9B2")
setEPS()
postscript("../latex/figures/grid_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:length(graph$loglike)), graph$loglike, type = "b", lty = 1, pch = 15, cex=.75, col = colors[1],
     xlab = "Iteration Number", ylab = "", ylim = c(-27, 10))
grid()
lines(c(1:length(graph$loglike)), graph$obj_fun, type = "b", xaxt = "n", lty = 2, pch=16, cex=.75, col = colors[2])
lines(c(1:length(graph$loglike)), graph$obj_fun - graph$loglike, type = "b", xaxt = "n", lty = 3, pch=17, cex=.75,
      col = colors[3])
legend("topright", legend = c("likelihood", "posterior", "prior"),
       col=colors, pch=c(15, 16, 17), lty=c(1, 2, 3), bty="n")
dev.off()
embed_fonts("../latex/figures/grid_trend.ps", outfile="../latex/figures/grid_trend.ps")
