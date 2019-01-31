library(igraph)
library(corrplot)
library(spectralGraphTopology)
library(pals)
library(extrafont)
set.seed(0)

n <- 64
p <- 30 * n
Strue <- wiener_kernel(n)
Ptrue <- solve(Strue)
Wtrue <- diag(diag(Ptrue)) - Ptrue
Y <- MASS::mvrnorm(p, mu = rep(0, n), Sigma = Strue)
graph <- learn_laplacian_matrix(cov(Y), k=1, w0 = "qp", beta = 1)
est_net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
true_net <- graph_from_adjacency_matrix(abs(Wtrue), mode = "undirected", weighted = TRUE)
blues <- brewer.blues(100)
reds <- brewer.reds(100)
blue_scale <- colorRamp(blues)
red_scale <- colorRamp(reds)
E(est_net)$color = apply(blue_scale(E(est_net)$weight / max(E(est_net)$weight)),
                         1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
non_zero_edges <- Ainv(Wtrue)[Ainv(Wtrue) != 0]
pos_edges <- non_zero_edges > 0
edge_colors <- blue_scale(pos_edges * abs(E(true_net)$weight) / max(pos_edges * abs(E(true_net)$weight)))
edge_colors <- 255 * edge_colors / max(edge_colors)
E(true_net)$color = apply(edge_colors, 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))

print(relativeError(graph$Lw, Ptrue))
print(graph$convergence)
print(graph$lambda)
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/wiener_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/wiener_true_graph.ps",
           family = "Times", height = 5, width = gr * 3.5)
plot(true_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/wiener_laplacian.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Lw / max(graph$Lw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
postscript("../latex/figures/wiener_covariance.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Strue / max(Strue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/wiener_precision.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Ptrue / max(Ptrue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()

