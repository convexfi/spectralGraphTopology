library(igraph)
library(viridis)
library(spectralGraphTopology)
library(glasso)

# design synthetic Laplacian of a grid graph
N <- 64
T <- 30 * N
grid <- make_lattice(length = sqrt(N), dim = 2)
# sample edges weights from U(.1, 3)
E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
Ltrue <- as.matrix(laplacian_matrix(grid))
# sample data from GP with covariance matrix set as
# the pseudo inverse of the true Laplacian
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
covY <- cov(Y)
# run spectralGraphTopology
graph <- learnGraphTopology(covY, beta = 1.8, Lwtol = 1e-4)
# run GLasso
graph_lasso <- glasso(covY, rho = 0.01)
# create an igraph object from the estimated and true adjacency matrices
est_net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
lasso_net <- graph_from_adjacency_matrix(diag(diag(graph_lasso$wi)) - graph_lasso$wi, mode = "undirected", weighted = TRUE)
true_net <- graph_from_adjacency_matrix(diag(diag(Ltrue)) - Ltrue, mode = "undirected", weighted = TRUE)
# colorify edges as per the estimated and true weights
colors <- inferno(5, begin = 0, end = 1, direction = -1)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(lasso_net)$color = apply(c_scale(abs(E(lasso_net)$weight) / max(abs(E(lasso_net)$weight))), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(true_net)$color = apply(c_scale(E(true_net)$weight / max(E(true_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
# create layout
la_est <- layout_on_grid(est_net)
lasso_est <- layout_on_grid(lasso_net)
la_true <- layout_on_grid(true_net)
# print numerical results
cat("\n** spectralGraphTopology results **\n")
cat("Relative error: ")
cat(relativeError(Ltrue, graph$Lw), "\n")
cat("Fscore: ")
cat(Fscore(Ltrue, graph$Lw, 1e-3), "\n")
cat("** glasso results **\n")
cat("Relative error: ")
cat(relativeError(Ltrue, graph_lasso$wi), "\n")
cat("Fscore: ")
cat(Fscore(Ltrue, graph_lasso$wi, 1e-3))
# plot graphs
plot(true_net, layout = la_true, vertex.label = NA, vertex.size = 3, main = "True Graph")
plot(est_net, layout = la_est, vertex.label = NA, vertex.size = 3, main = "spectralGraphTopology Graph")
plot(lasso_net, layout = lasso_est, vertex.label = NA, vertex.size = 3, main = "GLasso Graph")
