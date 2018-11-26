library(igraph)
library(spectralGraphTopology)
library(viridis)
library(corrplot)

set.seed(0)

N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
T <- as.integer(100 * N)
E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
Ltrue <- as.matrix(laplacian_matrix(grid))
Wtrue <- diag(diag(Ltrue)) - Ltrue
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
S <- cov(Y)

graph <- learnGraphTopology(S, w0 = "naive", beta = 4)
W <- diag(diag(graph$Lw)) - graph$Lw
Lnaive <- MASS::ginv(S)
Wnaive <- diag(diag(Lnaive)) - Lnaive
Wnaive[abs(Wnaive) < 1e-2] <- 0
grid_est <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
grid_naive <- graph_from_adjacency_matrix(Wnaive, mode = "undirected", weighted = TRUE)

colors <- inferno(5, begin = 0, end = 1, direction = -1)
c_scale <- colorRamp(colors)
E(grid_est)$color = apply(c_scale(E(grid_est)$weight / max(E(grid_est)$weight)), 1,
                          function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(grid)$color = apply(c_scale(E(grid)$weight / max(E(grid)$weight)), 1,
                      function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(grid_naive)$color = apply(c_scale(abs(E(grid_naive)$weight) / max(E(grid_naive)$weight)), 1,
                            function(x) rgb(x[1]/255, x[2]/255, x[3]/255))

la_est <- layout_on_grid(grid_est)
la_true <- layout_on_grid(grid)
la_naive <- layout_on_grid(grid_naive)

rel_spec <- relativeError(Ltrue, graph$Lw)
fs_spec <- Fscore(Ltrue, graph$Lw, 1e-2)
rel_naive <- relativeError(Ltrue, Lnaive)
fs_naive <- Fscore(Ltrue, Lnaive, 1e-2)
print(rel_spec)
print(fs_spec)
print(rel_naive)
print(fs_naive)

plot(grid, layout = la_true, vertex.label = NA, vertex.size = 3, main = "True Graph")
plot(grid_naive, layout = la_naive, vertex.label = NA, vertex.size = 3, main = "Naive Graph")
plot(grid_est, layout = la_est, vertex.label = NA, vertex.size = 3, main = "spectralGraphTopology Graph")
corrplot(Wtrue / max(Wtrue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
corrplot(W / max(W), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
corrplot(Wnaive / max(Wnaive), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
