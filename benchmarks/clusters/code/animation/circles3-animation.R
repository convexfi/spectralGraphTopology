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
S <- crossprod(t(circles3$data))
print(eigenvalues(S))
graph <- learn_k_component_graph(S, k = 3, beta = 1,
                                 maxiter = 570, record_weights = TRUE, record_objective = TRUE)
print(graph$obj_fun)
# pretty colors
colors <- c("#706FD3", "#FF5252", "#33D9B2")
gr = .5 * (1 + sqrt(5))
index <- c(1, 50, 100, 150, 200, 250, 300, 400, 500, length(graph$obj_fun))
j <- 0
for(i in index) {
  print(j)
  net <- graph_from_adjacency_matrix(A(c(as.matrix(graph$w_seq[[i]]))), mode = "undirected", weighted = TRUE)
  V(net)$cluster <- circles3$clusters
  E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                       function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                          colors[V(net)$cluster[x[1]]], 'grey'))
  V(net)$color <- colors[circles3$clusters]
  filename <- paste0("test", j, ".pdf")
  pdf(filename, height = 5, width = 3.5 * gr)
  plot(net, layout = circles3$data, vertex.label = NA, vertex.size = 3)
  dev.off()
  j <- j + 1
}

system("convert -dispose Background -delay 40 *.pdf circles3_reduced.gif")
file.remove(list.files(pattern=".pdf"))
