library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
set.seed(42)

eps <- 5e-2
N_realizations <- 1
ratios <- c(250)
# design synthetic Laplacian of a grid graph
N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
maxiter <- 5e4
for (j in 1:length(ratios)) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(grid))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_laplacian_matrix(S, w0 = "qp", beta = 1e3, fix_beta = TRUE,
                                    maxiter = maxiter, record_weights = TRUE,
                                    record_objective = TRUE)
    print(graph$convergence)
    niter <- length(graph$loglike)
    relative_error <- array(0, niter)
    fscore <- array(0, niter)
    for (i in 1:niter) {
      Lw <- L(as.array(graph$w_seq[[i]]))
      relative_error[i] <- relativeError(Ltrue, Lw)
      fscore[i] <- metrics(Ltrue, Lw, 5e-2)[1]
    }
  }
}

where_min_re <- c(1:niter)[relative_error == min(relative_error)]
# plot convergence trend
gr = .5 * (1 + sqrt(5))
colors <- c("#706FD3", "#FF5252", "#33D9B2")
pch <- c(15, 16, 17)
lty <- c(1, 2, 3)
setEPS()
postscript("grid_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
par(mar = c(5, 5, 3, 5))
plot(c(1:niter), graph$loglike, type = "b", pch = pch[1], cex=.3,
     ylim = c(min(graph$loglike, graph$obj_fun), max(graph$loglike, graph$obj_fun)),
     xlab = "Iteration Number", ylab = "Objective values", col = colors[1])
grid()
abline(v = where_min_re)
lines(c(1:niter), graph$obj_fun, type = "b", xaxt = "n", lty = 2, pch=16, cex=.3, col = colors[2])
par(new = TRUE)
plot(c(1:length(graph$loglike)), graph$obj_fun - graph$loglike, type = "b", pch = pch[3], cex=.3,
     xlab = "", ylab = "", xaxt = "n", lty = lty[3], yaxt = "n", col = colors[3])
axis(side = 4)
mtext("Constraint Value", side = 4, line = 3, family = "ComputerModern")
legend("right", legend = c("negloglike", "obj-func", "constr"), col=colors, pch = pch, lty = lty, bty="n")
dev.off()
embed_fonts("grid_trend.ps", outfile="grid_trend.ps")


setEPS()
postscript("relative_error.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:niter), relative_error, type = "b", lty = 1, pch = 15, cex=.2, col = "black",
     xlab = "Iteration Number", ylab = "Relative Error")
points(1, relative_error[1], type = "p", pch = 11, col = "red", cex = .7)
abline(v = where_min_re)
grid()
dev.off()
embed_fonts("relative_error.ps", outfile="relative_error.ps")

setEPS()
postscript("fscore.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:niter), fscore, type = "b", lty = 1, pch = 15, cex=.2, col = "black",
     xlab = "Iteration Number", ylab = "F - score")
points(1, fscore[1], type = "p", pch = 11, col = "red", cex = .7)
abline(v = where_min_re)
grid()
dev.off()
embed_fonts("fscore.ps", outfile="fscore.ps")
