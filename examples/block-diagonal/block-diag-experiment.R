library(igraph)
library(viridis)
library(corrplot)
library(spectralGraphTopology)
library(glasso)
set.seed(0)

N <- 64
T <- 30 * N
K <- 4
P <- diag(1, K)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
Ltrue <- as.matrix(laplacian_matrix(mgraph))
Wtrue <- diag(diag(Ltrue)) - Ltrue
# Erdo-Renyi as noise model
p <- .35
a <- .2
erdos_renyi <- erdos.renyi.game(N, p)
E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = 1)
Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
# Noisy Laplacian
Lnoisy <- (1 - a) * Ltrue + a * Lerdo
Wnoisy <- diag(diag(Lnoisy)) - Lnoisy
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Lnoisy))
graph <- learnGraphTopology(crossprod(Y) / T, K = 4, beta = 10)
est_net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
noisy_net <- graph_from_adjacency_matrix(Wnoisy, mode = "undirected", weighted = TRUE)
true_net <- graph_from_adjacency_matrix(Wtrue, mode = "undirected", weighted = TRUE)
colors <- blues9#viridis(5, begin = 0, end = 1, direction = -1)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(noisy_net)$color = apply(c_scale(E(noisy_net)$weight / max(E(noisy_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(true_net)$color = apply(c_scale(E(true_net)$weight / max(E(true_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
print(relativeError(Ltrue, graph$Lw))
print(Fscore(Ltrue, graph$Lw, 1e-3))
print(graph$convergence)
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("est_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("noisy_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(noisy_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("true_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(true_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("est_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$W / max(graph$W), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("noisy_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Wnoisy / max(Wnoisy), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("true_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Wtrue / max(Wtrue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
