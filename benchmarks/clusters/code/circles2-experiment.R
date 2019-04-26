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
S <- crossprod(t(circles2$data))
graph <- learn_k_component_graph(S, k = 2, beta = 1, fix_beta = FALSE, abstol = 1e-3)
# construct network
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
# use pretty colors for the nodes and edges
colors <- c("#706FD3", "#FF5252")
V(net)$cluster <- circles2$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#000000'))
V(net)$color <- colors[circles2$clusters]
# set up network plot
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/circles2.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, layout = circles2$data, vertex.label = NA, vertex.size = 3)
dev.off()
