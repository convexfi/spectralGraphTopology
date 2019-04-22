library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)
# number of points per cluster
N <- 100
# get datapoints
worms <- shapes.worms(N)
# learn underlying graph
graph <- learn_normalized_laplacian(cov(t(worms$data)), k = 2, maxiter = 1e4, beta = 1)
# build the network
Adjacency <- diag(diag(graph$NormalizedLaplacian)) - graph$NormalizedLaplacian
net <- graph_from_adjacency_matrix(Adjacency, mode = "undirected", weighted = TRUE)
# colorify edges and nodes
colors <- c("#706FD3", "#FF5252", "#33D9B2")
V(net)$cluster <- worms$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#000000'))
V(net)$color <- c(colors[1], colors[2])[worms$clusters]
# plot network
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/worms-normalized.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, layout = worms$data, vertex.label = NA, vertex.size = 3)
dev.off()
