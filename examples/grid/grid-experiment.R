library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(latex2exp)
library(colorspace)

N_realizations <- 10
ratios <- c(.5, .75, 1, 5, 10, 30, 100, 250, 500, 1000)
n_ratios <- c(1:length(ratios))
# design synthetic Laplacian of a grid graph
N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
rel_err_spec <- array(0, length(ratios))
rel_err_ggl <- array(0, length(ratios))
rel_err_gglA <- array(0, length(ratios))
rel_err_naive <- array(0, length(ratios))
fscore_spec <- array(0, length(ratios))
fscore_ggl <- array(0, length(ratios))
fscore_gglA <- array(0, length(ratios))
fscore_naive <- array(0, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab()
open(matlab)
print("success!")
A_mask <- matrix(1, 64, 64) - diag(64)
setVariable(matlab, A_mask = A_mask)
A_grid_mask <- as.matrix(laplacian_matrix(grid)) != 0
setVariable(matlab, A_grid_mask = A_grid_mask)

for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(grid))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- crossprod(Y) / T
    # run spectralGraphTopology
    if (ratios[j] <= 1)
      graph <- learnGraphTopology(S, w0 = "naive", beta = 1e-2, beta_max = 4, nbeta = 10)
    else
      graph <- learnGraphTopology(S, w0 = "naive", beta = 4)
    print(graph$convergence)

    # compute naive
    Lnaive <- MASS::ginv(S)
    # set data variable to MATLAB
    setVariable(matlab, S = S)
    evaluate(matlab, "[Lggl,~,~] = estimate_ggl(S, A_mask, 1e-3, 1e-4 , 1e-6, 40, 2)")
    evaluate(matlab, "[LgglA,~,~] = estimate_ggl(S, A_grid_mask, 1e-3, 1e-4 , 1e-6, 40, 2)")
    Lggl <- getVariable(matlab, "Lggl")
    LgglA <- getVariable(matlab, "LgglA")
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
    rel_spec = relativeError(Ltrue, graph$Lw)
    fs_spec = Fscore(Ltrue, graph$Lw, 1e-3)
    rel_ggl = relativeError(Ltrue, Lggl$Lggl)
    fs_ggl = Fscore(Ltrue, Lggl$Lggl, 1e-3)
    rel_gglA = relativeError(Ltrue, LgglA$LgglA)
    fs_gglA = Fscore(Ltrue, LgglA$LgglA, 1e-3)
    rel_naive = relativeError(Ltrue, Lnaive)
    fs_naive = Fscore(Ltrue, Lnaive, 1e-3)
    print(rel_spec)
    print(fs_spec)
    print(rel_ggl)
    print(fs_ggl)
    print(rel_gglA)
    print(fs_gglA)
    print(rel_naive)
    print(fs_naive)
    rel_err_spec[j] <- rel_err_spec[j] + rel_spec
    fscore_spec[j] <- fscore_spec[j] + fs_spec
    rel_err_ggl[j] <- rel_err_ggl[j] + rel_ggl
    fscore_ggl[j] <- fscore_ggl[j] + fs_ggl
    rel_err_gglA[j] <- rel_err_gglA[j] + rel_gglA
    fscore_gglA[j] <- fscore_gglA[j] + fs_gglA
    rel_err_naive[j] <- rel_err_naive[j] + rel_naive
    fscore_naive[j] <- fscore_naive[j] + fs_naive
  }
  rel_err_spec[j] <- rel_err_spec[j] / N_realizations
  fscore_spec[j] <- fscore_spec[j] / N_realizations
  rel_err_ggl[j] <- rel_err_ggl[j] / N_realizations
  fscore_ggl[j] <- fscore_ggl[j] / N_realizations
  rel_err_gglA[j] <- rel_err_gglA[j] / N_realizations
  fscore_gglA[j] <- fscore_gglA[j] / N_realizations
  rel_err_naive[j] <- rel_err_naive[j] / N_realizations
  fscore_naive[j] <- fscore_naive[j] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_spec[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_spec[j], "\n")
  cat("\n** GGL results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_ggl[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_ggl[j], "\n")
  cat("\n** GGL with A results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_gglA[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_gglA[j], "\n")
  cat("\n** Naive results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_naive[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_naive[j], "\n")
}

colors <- rainbow_hcl(4)
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error_grid.ps", family = "Times", height = 5, width = gr * 3.5)
plot(n_ratios, rel_err_naive, type = "b", pch=19, cex=.6,
     ylim=c(0, .6),
     xlab = TeX("$T / N$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
lines(n_ratios, rel_err_ggl, type = "b", pch=19, cex=.6, col = colors[2], xaxt = "n")
lines(n_ratios, rel_err_gglA, type = "b", pch=19, cex=.6, col = colors[3], xaxt = "n")
lines(n_ratios, rel_err_spec, type = "b", pch=19, cex=.6, col = colors[4], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("topright", legend = c(TeX("ISCM"), TeX("GGL"), TeX("GGL($\\mathbf{A}$)"), TeX("SGL")),
       col=colors, lty=c(1, 1, 1), cex=0.8)
dev.off()
setEPS()
postscript("fscore_grid.ps", family = "Times", height = 5, width = gr * 3.5)
plot(n_ratios, fscore_ggl, type = "b", pch=19, cex=.6,
     ylim=c(.5, 1.), xlab = TeX("$T / N$"), ylab = "Average F-score", col = colors[2], xaxt = "n")
lines(n_ratios, fscore_gglA, type = "b", pch=19, cex=.6, col = colors[3], xaxt = "n")
lines(n_ratios, fscore_spec, type = "b", pch=19, cex=.6, col = colors[4], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend = c(TeX("GGL"), TeX("GGL($\\mathbf{A}$)"), TeX("SGL")),
       col=colors[2:4], lty=c(1, 1, 1), cex=0.8)
dev.off()

#cat("** gggl results **\n")
#cat("Relative error: ")
#cat(relativeError(Ltrue, graph_ggl$wi), "\n")
#cat("Fscore: ")
#cat(Fscore(Ltrue, graph_ggl$wi, 1e-3))
# plot graphs
#plot(true_net, layout = la_true, vertex.label = NA, vertex.size = 3, main = "True Graph")
#plot(est_net, layout = la_est, vertex.label = NA, vertex.size = 3, main = "spectralGraphTopology Graph")
##plot(ggl_net, layout = ggl_est, vertex.label = NA, vertex.size = 3, main = "GLasso Graph")
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

