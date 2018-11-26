library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)
# number of points per cluster
N <- 100
circles2 <- shapes.circles2(N)
graph <- learnGraphTopology(crossprod(t(circles2$data)) + diag(rep(1/3, 2 * N)),
                            w0 = "naive", K = 2, beta = .2, ftol = 5e-5)
Lw <- graph$Lw
Lw <- diag(diag(Lw)) - Lw
net <- graph_from_adjacency_matrix(Lw, mode = "undirected", weighted = TRUE)
colors <- c("#706FD3", "#FF5252", "#33D9B2")
V(net)$cluster <- circles2$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#00000'))
V(net)$color <- c(colors[1], colors[2])[circles2$clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("circles2.ps", family = "Times", height = 5, width = gr * 3.5)
#plot(circles2$data, col=c(colors[1], colors[2])[circles2$clusters], pch = 19,
#     xlab = "X", ylab = "Y")
plot(net, layout = circles2$data, vertex.label = NA, vertex.size = 3)
dev.off()
