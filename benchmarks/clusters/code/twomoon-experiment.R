library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)
# number of nodes per cluster
N <- 100
# generate datapoints
twomoon <- shapes.two.moon(N)
# estimate underlying graph
S <- crossprod(t(twomoon$data))
graph <- learn_laplacian_matrix(S, k = 2, beta = 1, tol = 1e-2)
c(graph$eigenvalues)
graph$beta_seq
# build network
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
# colorify nodes and edges
colors <- c("#706FD3", "#FF5252", "#33D9B2")
V(net)$cluster <- twomoon$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                      function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#000000'))
V(net)$color <- c(colors[1], colors[2])[twomoon$clusters]
# plot network
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/twomoon.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, layout = twomoon$data, vertex.label = NA, vertex.size = 3)
dev.off()
