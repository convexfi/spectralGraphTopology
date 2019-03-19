library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

eps <- 1e-2
n_realizations <- 500
ratios <- c(10)
n <- 16
k <- 4
P <- diag(1, k)
mgraph <- sample_sbm(n, pref.matrix = P, block.sizes = c(rep(n / k, k)))
maxiter <- 5e4
relative_error_list <- list()
fscore_list <- list()
nll_list <- list()
objfun_list <- list()
constr_list <- list()

for (j in 1:length(ratios)) {
  t <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", t, "samples per node, t/n = ", ratios[j], "\n")
  for (r in 1:n_realizations) {
    print(r)
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(t, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_laplacian_matrix(S, w0 = "naive", k = 4, beta = 1e2, fix_beta = TRUE,
                                    edge_tol = eps, maxiter = maxiter,
                                    record_weights = TRUE, record_objective = TRUE)
    niter <- length(graph$loglike)
    relative_error <- array(0, niter)
    fscore <- array(0, niter)
    for (i in 1:niter) {
      Lw <- L(as.array(graph$w_seq[[i]]))
      relative_error[i] <- relativeError(Ltrue, Lw)
      fscore[i] <- metrics(Ltrue, Lw, eps)[1]
    }
    relative_error_list <- rlist::list.append(relative_error_list, relative_error)
    fscore_list <- rlist::list.append(fscore_list, fscore)
    nll_list <- rlist::list.append(nll_list, graph$loglike)
    objfun_list <- rlist::list.append(objfun_list, graph$obj_fun)
    constr_list <- rlist::list.append(constr_list, graph$obj_fun - graph$loglike)
  }
}

saveRDS(relative_error_list, file = "relerr_beta1e2.RDS")
saveRDS(fscore_list, file = "fscore_beta1e2.RDS")
saveRDS(nll_list, file = "nll_beta1e2.RDS")
saveRDS(objfun_list, file = "objfun_beta1e2.RDS")
saveRDS(constr_list, file = "constr_beta1e2.RDS")

# plot convergence trend
#gr = .5 * (1 + sqrt(5))
#colors <- c("#706FD3", "#FF5252", "#33D9B2")
#pch <- c(15, 16, 17)
#lty <- c(1, 2, 3)
#setEPS()
#postscript("block_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#par(mar = c(5, 5, 3, 5))
#plot(c(1:niter), graph$loglike, type = "l",
#     ylim = c(min(graph$loglike, graph$obj_fun), -4), lty = lty[1],
#     xlab = "Iteration Number", ylab = "Objective values", col = colors[1])
#lines(c(1:niter), graph$obj_fun, type = "l", xaxt = "n", lty = lty[2], col = colors[2])
#par(new = TRUE)
#plot(c(1:length(graph$loglike)), graph$obj_fun - graph$loglike, ylim = c(0, 10), type = "l",
#     xlab = "", ylab = "", xaxt = "n", lty = lty[3], yaxt = "n", col = colors[3])
#axis(side = 4)
#mtext("Constraint Value", side = 4, line = 3, family = "ComputerModern")
#legend("topright", legend = c("Negative Loglikelihood", "Objective Function", "Constraint Value"),
#       col=colors, lty = lty, bty="n")
#dev.off()
#embed_fonts("block_trend.ps", outfile="block_trend.ps")
#
#setEPS()
#postscript("relative_error.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(c(1:niter), relative_error, type = "l", lty = 1, col = "black",
#     xlab = "Iteration Number", ylab = "Relative Error")
#points(1, relative_error[1], pch = 11, cex = .5, col = "red")
#dev.off()
#embed_fonts("relative_error.ps", outfile="relative_error.ps")
#
#setEPS()
#postscript("fscore.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(c(1:niter), fscore, type = "l", lty = 1, col = "black",
#     xlab = "Iteration Number", ylab = "F - score")
#points(1, fscore[1], pch = 11, cex = .5, col = "red")
#dev.off()
#embed_fonts("fscore.ps", outfile="fscore.ps")
