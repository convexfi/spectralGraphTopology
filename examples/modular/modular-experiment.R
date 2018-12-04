library(igraph)
library(viridis)
library(spectralGraphTopology)

set.seed(0)
N <- 64
#T <- 30 * N
T <- 100 * N
P <- diag(.49, 4) + 0.01
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(16, 4)))
E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 1)
Ltrue <- as.matrix(laplacian_matrix(mgraph))
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
graph <- learnGraphTopology(cov(Y), K = 1, w0 = "qp", beta = 5)
W <- graph$W
W[W < 1e-1] = 0
est_net <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
true_net <- graph_from_adjacency_matrix(diag(diag(Ltrue)) - Ltrue, mode = "undirected", weighted = TRUE)
colors <- inferno(5, begin = 0, end = 1, direction = -1)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(true_net)$color = apply(c_scale(E(true_net)$weight / max(E(true_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
print(relativeError(Ltrue, graph$Lw))
print(Fscore(Ltrue, graph$Lw, 1e-1))
print(graph$convergence)
plot(est_net, vertex.label = NA, vertex.size = 3)
plot(true_net, vertex.label = NA, vertex.size = 3)
plot(graph$elapsed_time, graph$obj_fun, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
     xlab = "CPU time [seconds]", ylab = "Objective function")
