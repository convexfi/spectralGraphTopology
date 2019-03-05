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
graph <- learn_laplacian_matrix(crossprod(t(circles3$data)), k = 3, beta = 1, tol = 1e-2, record_weights = TRUE)
print(length(graph$beta_seq))
# pretty colors
colors <- c("#706FD3", "#FF5252", "#33D9B2")
# colorify edges and nodes
# construct the network
# plot
gr = .5 * (1 + sqrt(5))
pdf("circles3.pdf", height = 5, width = 3.5 * gr)
for(i in c(length(graph$beta_seq))) {
  net <- graph_from_adjacency_matrix(A(c(as.matrix(graph$w_seq[[i]]))), mode = "undirected", weighted = TRUE)
  V(net)$cluster <- circles3$clusters
  E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                       function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                          colors[V(net)$cluster[x[1]]], '#000000'))
  V(net)$color <- colors[circles3$clusters]
  plot(net, layout = circles3$data, vertex.label = NA, vertex.size = 3)
}
dev.off()

system("convert -delay 40 *.pdf circles3_reduced.gif")
file.remove(list.files(pattern=".pdf"))
