library(igraph)
library(viridis)
library(spectralGraphTopology)

N_realizations <- 100
ratios <- c(.5, .75, 1, 5, 10, 30, 100, 250, 500, 1000)
#ratios <- c(100, 250, 500, 1000)
# design synthetic Laplacian of a grid graph
N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
# sample edges weights from U(.1, 3)
rel_err <- array(0, length(ratios))
fscore <- array(0, length(ratios))
for (j in c(1:length(ratios))) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(grid))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    # run spectralGraphTopology
    if (ratios[j] <= 1)
      graph <- learnGraphTopology(crossprod(Y) / T, w0 = "naive",
                                  beta = 1e-2, beta_max = 4, npoints = 5)
    else
      graph <- learnGraphTopology(crossprod(Y) / T, w0 = "naive", beta = 4)
    print(graph$convergence)
    # run GLasso
    # graph_lasso <- glasso(crossprod(Y) / T, rho = 5e-3)
    # create an igraph object from the estimated and true adjacency matrices
    #est_net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
    ##lasso_net <- graph_from_adjacency_matrix(diag(diag(graph_lasso$wi)) - graph_lasso$wi, mode = "undirected", weighted = TRUE)
    #true_net <- graph_from_adjacency_matrix(diag(diag(Ltrue)) - Ltrue, mode = "undirected", weighted = TRUE)
    ## colorify edges as per the estimated and true weights
    #colors <- inferno(5, begin = 0, end = 1, direction = -1)
    #c_scale <- colorRamp(colors)
    #E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
    ##E(lasso_net)$color = apply(c_scale(abs(E(lasso_net)$weight) / max(abs(E(lasso_net)$weight))), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
    #E(true_net)$color = apply(c_scale(E(true_net)$weight / max(E(true_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
    ## create layout
    #la_est <- layout_on_grid(est_net)
    ##lasso_est <- layout_on_grid(lasso_net)
    #la_true <- layout_on_grid(true_net)
    # print numerical results
    rel = relativeError(Ltrue, graph$Lw)
    fs = Fscore(Ltrue, graph$Lw, 1e-3)
    print(rel)
    print(fs)
    #print(relativeError(Ltrue, graph_lasso$wi))
    rel_err[j] <- rel_err[j] + rel
    fscore[j] <- fscore[j] + fs
  }
  rel_err[j] <- rel_err[j] / N_realizations
  fscore[j] <- fscore[j] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore[j], "\n")
}

gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error_grid.ps", family = "Times", height = 5, width = gr * 3.5)
plot(c(1:length(ratios)), rel_err, type = "b", pch=19, cex=.6,
     ylim=c(min(rel_err) - 1e-2, max(rel_err) + .1e-2),
     xlab = "T / N", ylab = "Average Relative Error", col = "black", xaxt = "n")
axis(side = 1, at = c(1:length(ratios)), labels = ratios)
dev.off()
setEPS()
postscript("fscore_grid.ps", family = "Times", height = 5, width = gr * 3.5)
plot(c(1:length(ratios)), fscore, type = "b", pch=19, cex=.6,
     ylim=c(min(fscore) - 1e-2, 1.),
     xlab = "T / N", ylab = "Average F-score", col = "black", xaxt = "n")
axis(side = 1, at = c(1:length(ratios)), labels = ratios)
dev.off()

#cat("** glasso results **\n")
#cat("Relative error: ")
#cat(relativeError(Ltrue, graph_lasso$wi), "\n")
#cat("Fscore: ")
#cat(Fscore(Ltrue, graph_lasso$wi, 1e-3))
# plot graphs
#plot(true_net, layout = la_true, vertex.label = NA, vertex.size = 3, main = "True Graph")
#plot(est_net, layout = la_est, vertex.label = NA, vertex.size = 3, main = "spectralGraphTopology Graph")
##plot(lasso_net, layout = lasso_est, vertex.label = NA, vertex.size = 3, main = "GLasso Graph")
#plot(graph$elapsed_time, graph$obj_fun, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
#     xlab = "CPU time [seconds]", ylab = "Objective function")

#relerr_seq <- c()
#ii <- c(1:length(graph$w_seq))
#for (i in ii) {
#  relerr <- relativeError(Ltrue, L(graph$w_seq[[i]]))
#  relerr_seq <- c(relerr_seq, relerr)
#}
#plot(graph$elapsed_time, relerr_seq, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
#     xlab = "CPU time [seconds]", ylab = "Relative Error")

