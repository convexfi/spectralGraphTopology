library(kernlab)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)
data(spirals)

# learn underlying graph
n <- nrow(spirals)
graph <- learn_k_component_graph(crossprod(t(spirals)), k = 2, beta = 1, tol = 1e-2)
# construct network
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
clusters <- components(net)$membership
colors <- c("#706FD3", "#FF5252")
V(net)$cluster <- clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#000000'))
V(net)$color <- colors[clusters]
# use pretty colors for the nodes and edges
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/spirals.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, layout = spirals, vertex.label = NA, vertex.size = 3)
dev.off()
