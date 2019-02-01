library(spectralGraphTopology)
library(corrplot)
library(pals)
library(extrafont)
library(igraph)
set.seed(42)
n1 <- 8
n2 <- 4
bipartite <- sample_bipartite(n1, n2, type="Gnp", p = .9, directed=FALSE)
E(bipartite)$weight <- runif(gsize(bipartite), min = 1, max = 2)

#erdos_renyi <- erdos.renyi.game(n1 + n2, .3)
#E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = .5)
#Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))

Lw <- as.matrix(laplacian_matrix(bipartite))
Aw <- diag(diag(Lw)) - Lw
w_true <- Linv(Lw)
n <- ncol(Lw)
p <- 10 * n
# get datapoints
Y <- MASS::mvrnorm(p, rep(0, n), Sigma = MASS::ginv(Lw))
# learn underlying graph
S <- cov(Y)
graph <- learn_adjacency_and_laplacian(S, w0 = "naive", beta1 = 10, beta2 = .5, maxiter = 1e5)
print(graph$lambda)
print(graph$psi)
print(graph$obj_fun)
Sinv <- MASS::ginv(S)
w_naive <- spectralGraphTopology:::w_init(w0 = "naive", Sinv)
w_qp <- spectralGraphTopology:::w_init(w0 = "qp", Sinv)
print(relativeError(graph$Aw, Aw))
print(Fscore(graph$Aw, Aw, 5e-1))
print(relativeError(A(w_qp), Aw))
print(Fscore(A(w_qp), Aw, 5e-1))
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
legend("topright", legend = c("objective-function"),
       col=c("#706FD3"), pch=c(15), lty=c(1), bty="n")
dev.off()
embed_fonts("test_trend.ps", outfile="test_trend.ps")
