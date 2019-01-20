library(spectralGraphTopology)
library(corrplot)
library(pals)
library(extrafont)
library(igraph)
set.seed(42)

bipartite <- sample_bipartite(24, 8, type="Gnp", p = .4, directed=FALSE)
E(bipartite)$weight <- runif(gsize(bipartite), min = 1e-1, max = 1)
Lw <- as.matrix(laplacian_matrix(bipartite))
Aw <- diag(diag(Lw)) - Lw
w_true <- Linv(Lw)
n <- ncol(Lw)
T <- 100 * n
# get datapoints
Y <- MASS::mvrnorm(T, rep(0, n), Sigma = MASS::ginv(Lw))
# learn underlying graph
graph <- learnAdjacencyGraphTopology(cov(Y), w0 = "qp", beta = 1, beta_max = 100, nbeta = 20, maxiter = 1e5)
print(graph)
print(relativeError(graph$Aw, Aw))
print(Fscore(graph$Aw, Aw, 5e-2))
# build the network
net <- graph_from_adjacency_matrix(Aw, mode = "undirected", weighted = TRUE)
V(net)$type = V(bipartite)$type
# plot network
colors <- brewer.blues(20)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(bipartite)$color = apply(c_scale(E(bipartite)$weight / max(E(bipartite)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
gr = .5 * (1 + sqrt(5))
#setEPS()
#postscript("test.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(net, vertex.label = NA,
#     vertex.color=c("#706FD3", "#33D9B2")[V(net)$type+1],
#     vertex.size = 3)
#dev.off()
#embed_fonts("test.ps", outfile="test.ps")
#
#setEPS()
#postscript("test_true.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(bipartite, vertex.label = NA,
#     vertex.color=c("#706FD3", "#33D9B2")[V(bipartite)$type+1],
#     vertex.size = 3)
#dev.off()
#embed_fonts("test_true.ps", outfile="test_true.ps")

setEPS()
postscript("est_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Aw / max(graph$Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("true_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Aw / max(Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()

# plot convergence trend
niter <- length(graph$loglike)
iterations <- c(1:niter)
setEPS()
postscript("test_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(iterations, graph$obj_fun, type = "b", lty = 1, pch = 15, cex=.75, col = "#706FD3",
     xlab = "Iteration number", ylab = "Objective value")
grid()
lines(iterations, graph$loglike, type = "b", xaxt = "n", lty = 2, pch=16, cex=.75, col = "#33D9B2")
legend("topright", legend = c("posterior", "loglikelihood"),
       col=c("#706FD3", "#33D9B2"), pch=c(15, 16), lty=c(1, 2), bty="n")
dev.off()
embed_fonts("test_trend.ps", outfile="test_trend.ps")
