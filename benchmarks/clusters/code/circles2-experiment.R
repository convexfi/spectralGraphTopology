library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)
# number of points per cluster
N <- 100
# generate nodes
circles2 <- shapes.circles2(N)
# learn underlying graph
graph <- learn_laplacian_matrix(crossprod(t(circles2$data)) + diag(rep(1/3, 2 * N)),
                                w0 = "naive", k = 2, beta = .25, ftol = 1e-4, Lwtol = 1e-4)
# construct network
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
# use pretty colors for the nodes and edges
colors <- c("#706FD3", "#FF5252", "#33D9B2")
V(net)$cluster <- circles2$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#000000'))
V(net)$color <- c(colors[1], colors[2])[circles2$clusters]
# set up network plot
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/circles2.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, layout = circles2$data, vertex.label = NA, vertex.size = 3)
dev.off()
