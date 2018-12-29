library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)
# number of points per cluster
N <- 100
# generate datapoints
circles3 <- shapes.circles3(N)
# learn underlying graph
graph <- learnGraphTopology(crossprod(t(circles3$data)) + diag(rep(1/3, 3 * N)),
                            w0 = "naive", K = 3, beta = .25, ftol = 1e-3, Lwtol = 1e-3)
# construct the network
net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
# pretty colors
colors <- c("#706FD3", "#FF5252", "#33D9B2")
# colorify edges and nodes
V(net)$cluster <- circles3$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#00000'))
V(net)$color <- c(colors[1], colors[2], colors[3])[circles3$clusters]
# plot
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/circles3.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, layout = circles3$data, vertex.label = NA, vertex.size = 3)
dev.off()
