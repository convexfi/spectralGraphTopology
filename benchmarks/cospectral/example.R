library(spectralGraphTopology)
library(corrplot)
library(igraph)
library(latex2exp)
library(pals)
set.seed(123)

p <- 15
f <- .05
w <- runif(p * (p - 1) / 2)
Laplacian <- L(w)
lambda <- spectralGraphTopology:::eigval_sym(Laplacian)
Adjacency <- diag(diag(Laplacian)) - Laplacian
# construct the network
true_net <- graph_from_adjacency_matrix(Adjacency, mode = "undirected", weighted = TRUE)
# pretty colors
colors <- c("#2d3436")
colors_edge <- brewer.blues(20)
c_scale <- colorRamp(colors_edge)
# colorify edges and nodes
V(true_net)$cluster <- c(rep(1, p))
E(true_net)$color = apply(c_scale(E(true_net)$weight / max(E(true_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(true_net)$color <- colors[c(rep(1, p))]

# plot
l <- layout.random(true_net)
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("original-graph.ps", height = 5, fonts = "Palatino", width = gr * 3.5)
plot(true_net, vertex.label = NA, vertex.size = 5, vertex.color = "#686de0", layout = l)
dev.off()

Y <- MASS::mvrnorm(p * 50, rep(0, p), MASS::ginv(Laplacian))
S <- cov(Y)
graph <- learn_cospectral_graph(S, lambda = f * lambda[2:p], beta = 100)
graph_est <- learn_k_component_graph(S, beta = 100)

graph$Adjacency[graph$Adjacency < 3e-2] <- 0
graph_est$Adjacency[graph_est$Adjacency < 3e-2] <- 0
print(sum(spectralGraphTopology:::Ainv(Adjacency) > 0))
print(sum(spectralGraphTopology:::Ainv(graph$Adjacency) > 0))
print(sum(spectralGraphTopology:::Ainv(graph_est$Adjacency) > 0))
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
V(net)$cluster <- c(rep(1, p))
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color <- colors[c(rep(1, p))]
setEPS()
postscript("estimated-graph-lower-lambda.ps", height = 5, fonts = "Palatino", width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 5, vertex.color = "#686de0", layout = l)
dev.off()

# estimated original graph
net <- graph_from_adjacency_matrix(graph_est$Adjacency, mode = "undirected", weighted = TRUE)
V(net)$cluster <- c(rep(1, p))
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color <- colors[c(rep(1, p))]
setEPS()
postscript("estimated-graph-original-lambda.ps", height = 5, fonts = "Palatino", width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 5, vertex.color = "#686de0", layout = l)
dev.off()

eigvals_label_fun <- function(k) {
  labels <- c()
  for (i in 1:k)
    labels <- c(labels, TeX(paste0("$\\lambda_{", toString(i), "}$")))
  return(labels)
}

eigvals_labels <- eigvals_label_fun(p)
spectralGraphTopology:::metrics(graph_est$Laplacian, Laplacian, 1e-3)
relative_error(graph_est$Laplacian, Laplacian)

setEPS()
postscript("eigenvalues.ps", height = 5, fonts = "Palatino", width = gr * 3.5)
par(family = "Palatino")
plot(c(1:p), f * lambda,
     ylab = "", xlab = "Eigenvalues", xaxt = "n", pch = 17, col = "#ff7675")
points(c(1:p), abs(spectralGraphTopology:::eigval_sym(graph$Laplacian)),
       pch = 11, col = "#0984e3")
axis(side = 1, at = c(1:p), labels = eigvals_labels)
grid()
legend("topleft", legend=c("Ground-truth", "Estimated"),
       col=c("#ff7675", "#0984e3"), pch=c(17, 11))
dev.off()
