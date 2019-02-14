library(spectralGraphTopology)
library(corrplot)
library(pals)
library(extrafont)
library(igraph)
library(clusterSim)
library(scales)
set.seed(234)

n1 <- 40
n2 <- 24
n <- n1 + n2
pc <- .7
p <- 500 * n

bipartite <- sample_bipartite(n1, n2, type="Gnp", p = pc, directed=FALSE)
# randomly assign edge weights to connected nodes
E(bipartite)$weight <- runif(gsize(bipartite), min = .1, max = 1)
# Erdo-Renyi as noise model
a <- .5
erdos_renyi <- erdos.renyi.game(n, p = .35)
E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
# get true Laplacian and Adjacency
Lw <- as.matrix(laplacian_matrix(bipartite))
Aw <- diag(diag(Lw)) - Lw
w_true <- Linv(Lw)
Lnoisy = Lw + Lerdo
Anoisy <- diag(diag(Lnoisy)) - Lnoisy
# set number of samples
Y <- MASS::mvrnorm(p, rep(0, n), Sigma = MASS::ginv(Lnoisy))
S <- cov(Y)
graph <- learn_bipartite_graph(S, z = abs(n2 - n1), w0 = "qp", beta = 1e5, ftol = 1e-6, maxiter = 1e5)
#Sinv <- MASS::ginv(S)
#w_naive <- spectralGraphTopology:::w_init(w0 = "naive", Sinv)
#w_qp <- spectralGraphTopology:::w_init(w0 = "qp", Sinv)
#Lnaive <- L(w_naive)
#Lqp <- L(w_qp)
#Aqp <- diag(diag(Lqp)) - Lqp
print(relativeError(Aw, graph$Aw))
print(Fscore(Aw, graph$Aw, 1e-2))
#print(relativeError(Aw, Lqp))
#print(Fscore(Aw, Aqp, 1e-2))

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
postscript("est_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Aw / max(graph$Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("noisy_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Anoisy / max(Anoisy), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("true_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Aw / max(Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()

## build the network
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
#V(net)$type = V(bipartite)$type
## plot network
colors <- alpha(brewer.greys(100), alpha = .1)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
clusters <- c(rep(1, n1), rep(2, n2))
V(net)$color <- c("#706FD3", "#33D9B2")[clusters]
V(net)$shape <- c("circle", "square")[clusters]
setEPS()
postscript("test.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
y <- c(rep(.2, n1), rep(.9, n2))
x <- c(.5 * c(1:n1) - .4, .5 * c(1:n2))
plot(net, vertex.label = NA, layout = t(rbind(y, x)), vertex.size = 4)
dev.off()
embed_fonts("test.ps", outfile="test.ps")
