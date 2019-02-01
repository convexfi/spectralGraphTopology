library(igraph)
library(corrplot)
library(spectralGraphTopology)
library(pals)
library(huge)
library(extrafont)
set.seed(0)

n <- 32
p <- 100*n
Strue <- periodic_kernel(2, sqrt(2), 10, n) + gaussian_kernel(2, 5e-1, n)
Ptrue <- spectralGraphTopology:::inv_sympd(Strue)
Wtrue <- diag(diag(Ptrue)) - Ptrue
true_net <- graph_from_adjacency_matrix(abs(Wtrue), mode = "undirected", weighted = TRUE)
Y <- MASS::mvrnorm(p, mu = rep(0, n), Sigma = Strue)
graph <- learn_laplacian_matrix(cov(Y), k = 1, w0 = "qp", beta = 20, maxiter = 1e5)
glasso <- huge(cov(Y), method = "glasso")
Pglasso <- matrix(glasso$icov[[10]], nrow = 32)
est_net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
blues <- brewer.blues(100)
reds <- brewer.reds(100)
blue_scale <- colorRamp(blues)
red_scale <- colorRamp(reds)
E(est_net)$color = apply(blue_scale(E(est_net)$weight / max(E(est_net)$weight)),
                         1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
non_zero_edges <- Ainv(Wtrue)[Ainv(Wtrue) != 0]
pos_edges <- non_zero_edges > 0
neg_edges <- !pos_edges
edge_colors <- (blue_scale(pos_edges * E(true_net)$weight / max(pos_edges * E(true_net)$weight)) +
                red_scale(neg_edges * abs(E(true_net)$weight) / max(neg_edges * abs(E(true_net)$weight))))
edge_colors <- 255 * edge_colors / max(edge_colors)
E(true_net)$color = apply(edge_colors, 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
print(graph$convergence)
print(graph$lambda)

gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/periodic_gaussian_graph.ps",
           family = "Times", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_true_graph.ps",
           family = "Times", height = 5, width = gr * 3.5)
plot(true_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_laplacian.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Lw / max(graph$Lw), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_glasso_precision.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(Pglasso / max(Pglasso), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
postscript("../latex/figures/periodic_gaussian_covariance.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(Strue / max(Strue), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_precision.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(Ptrue / max(Ptrue), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
