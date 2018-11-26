library(igraph)
library(viridis)
library(corrplot)
library(spectralGraphTopology)
library(pals)
set.seed(0)

N <- 12
T <- 10 * N
K <- 4
P <- diag(1, K)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
Ltrue <- as.matrix(laplacian_matrix(mgraph))
Wtrue <- diag(diag(Ltrue)) - Ltrue
# Erdo-Renyi as noise model
p <- 1.
a <- .1
erdos_renyi <- erdos.renyi.game(N, p)
E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
# Noisy Laplacian
Lnoisy <- Ltrue + Lerdo
Wnoisy <- diag(diag(Lnoisy)) - Lnoisy
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Lnoisy))
graph <- learnGraphTopology(crossprod(Y) / T, K, w0 = "naive", beta = .5)
ULUT <- crossprod(sqrt(graph$lambda) * t(graph$U))
W <- abs(diag(diag(ULUT)) - ULUT)
est_net <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
noisy_net <- graph_from_adjacency_matrix(Wnoisy, mode = "undirected", weighted = TRUE)
true_net <- graph_from_adjacency_matrix(Wtrue, mode = "undirected", weighted = TRUE)
colors <- brewer.greys(20)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(noisy_net)$color = apply(c_scale(E(noisy_net)$weight / max(E(noisy_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(true_net)$color = apply(c_scale(E(true_net)$weight / max(E(true_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))

print(relativeError(Ltrue, ULUT))
print(Fscore(Ltrue, ULUT, 1e-3))
print(graph$convergence)
print(graph$beta)

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
corrplot(W / max(W), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("noisy_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Wnoisy / max(Wnoisy), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("true_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Wtrue / max(Wtrue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
