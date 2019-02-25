library(kernlab)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)
data(spirals)

# learn underlying graph
n <- nrow(spirals)
graph <- learn_laplacian_matrix(crossprod(t(spirals)) + diag(1/3, n), w0 = "naive", k = 2,
                                beta = .15, alpha = 1e-2, ftol = 1e-4, Lwtol = 1e-4)
# construct network
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
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
