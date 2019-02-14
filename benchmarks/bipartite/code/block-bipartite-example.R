library(spectralGraphTopology)
library(corrplot)
library(pals)
library(extrafont)
library(igraph)
library(clusterSim)
library(scales)
set.seed(234)

n1 <- 10
n2 <- 4
n <- 3 * n1 + 3 * n2 - 10
p <- 250 * n

# randomly assign edge weights to connected nodes
b1 <- sample_bipartite(n1, n2, type="Gnp", p = .7, directed=FALSE)
b2 <- sample_bipartite(n1-4, n2, type="Gnp", p = .8, directed=FALSE)
b3 <- sample_bipartite(n1-6, n2, type="Gnp", p = .9, directed=FALSE)
E(b1)$weight <- runif(gsize(b1), min = 1, max = 3)
E(b2)$weight <- runif(gsize(b2), min = 1, max = 3)
E(b3)$weight <- runif(gsize(b3), min = 1, max = 3)
Lw1 <- as.matrix(laplacian_matrix(b1))
Aw1 <- diag(diag(Lw1)) - Lw1
Lw2 <- as.matrix(laplacian_matrix(b2))
Aw2 <- diag(diag(Lw2)) - Lw2
Lw3 <- as.matrix(laplacian_matrix(b3))
Aw3 <- diag(diag(Lw3)) - Lw3
Lw <- blockDiag(Lw1, Lw2, Lw3)
Aw <- blockDiag(Aw1, Aw2, Aw3)

# Erdo-Renyi as noise model
a <- 1
erdos_renyi <- erdos.renyi.game(n, p = .35)
E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
Lnoisy <- Lw + Lerdo
Anoisy <- diag(diag(Lnoisy)) - Lnoisy
# set number of samples
Y <- MASS::mvrnorm(p, rep(0, n), Sigma = MASS::ginv(Lnoisy))
S <- cov(Y)
graph <- learn_adjacency_and_laplacian(S, k = 3, z = 0, w0 = "qp", beta1 = 1e4,
                                       beta2 = 1e5, alpha = 1e-2)
print(graph$convergence)
print(relativeError(Aw, graph$Aw))
print(Fscore(Aw, graph$Aw, 1e-1))

## build the network
#net <- graph_from_adjacency_matrix(Aw, mode = "undirected", weighted = TRUE)
#V(net)$type = V(bipartite)$type
## plot network
#colors <- brewer.blues(20)
#c_scale <- colorRamp(colors)
#E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
#E(bipartite)$color = apply(c_scale(E(bipartite)$weight / max(E(bipartite)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
gr = .5 * (1 + sqrt(5))
#setEPS()
#postscript("test.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(net, vertex.label = NA, layout = layout_as_bipartite,
#     vertex.color=c("#706FD3", "#33D9B2")[V(net)$type+1],
#     vertex.size = 3)
#dev.off()
#embed_fonts("test.ps", outfile="test.ps")
#
#setEPS()
#postscript("test_true.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(bipartite, vertex.label = NA, layout = layout_as_bipartite,
#     vertex.color=c("#706FD3", "#33D9B2")[V(bipartite)$type+1],
#     vertex.size = 3)
#dev.off()
#embed_fonts("test_true.ps", outfile="test_true.ps")
#
setEPS()
postscript("../latex/figures/est_block_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Aw / max(graph$Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/true_block_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Aw / max(Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/noisy_block_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Anoisy / max(Anoisy), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()


## build the network
graph$Aw[graph$Aw < 1e-1] <- 0
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
#V(net)$type = V(bipartite)$type
## plot network
colors <- c("#706FD3", "#33D9B2")
clusters <- c(rep(1, n1), rep(2, n2), rep(1, n1-4), rep(2, n2), rep(1, n1-6), rep(2, n2))
V(net)$cluster <- clusters
V(net)$color <- colors[clusters]
V(net)$shape <- c("circle", "square")[clusters]
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                      function(x) ifelse(V(net)$cluster[x[1]] != V(net)$cluster[x[2]],
                                         "pink", brewer.greys(5)[2]))
setEPS()
postscript("../latex/figures/est_block_graph.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
y <- c(rep(.2, n1), rep(.5, n2), rep(.2, n1-4), rep(.5, n2), rep(.2, n1-6), rep(.5, n2))
x <- c(1:(3 * (n1 + n2) - 10))
plot(net, vertex.label = NA, layout = t(rbind(x, y)), vertex.size = 4)
#plot(net, vertex.label = NA, vertex.size = 4)
dev.off()
embed_fonts("../latex/figures/est_block_graph.ps", outfile="../latex/figures/est_block_graph.ps")
## plot convergence trend
#niter <- length(graph$loglike)
#iterations <- c(1:niter)
#setEPS()
#postscript("test_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(iterations, graph$obj_fun, type = "b", lty = 1, pch = 15, cex=.75, col = "#706FD3",
#     xlab = "Iteration number", ylab = "Objective value")
#grid()
#legend("topright", legend = c("objective-function"),
#       col=c("#706FD3"), pch=c(15), lty=c(1), bty="n")
#dev.off()
#embed_fonts("test_trend.ps", outfile="test_trend.ps")
