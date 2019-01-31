library(igraph)
library(corrplot)
library(spectralGraphTopology)
library(pals)
library(extrafont)
set.seed(0)

n <- 512
p <- 100 * n
k <- 8
WK <- wiener_kernel(n/8)
Strue <- blockDiag(WK, WK, WK, WK,
                   WK, WK, WK, WK)
Y <- MASS::mvrnorm(p, mu = rep(0, n), Sigma = Strue)
graph <- learn_laplacian_matrix(cov(Y), k=k, w0 = "naive", beta = 4)
est_net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
colors <- brewer.blues(20)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))

print(graph$convergence)
print(graph$lambda)

gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/wiener_block_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/wiener_block_laplacian.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Lw / max(graph$Lw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
postscript("../latex/figures/wiener_block_covariance.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Strue / max(Strue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/wiener_block_precision.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(solve(Strue) / max(solve(Strue)), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()

