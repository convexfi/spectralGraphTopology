library(corrplot)
library(spectralGraphTopology)
library(igraph)
library(pals)

S <- readRDS("correlation-matrix.rds")
graph <- learn_laplacian_matrix(S, k = 40, beta = 1e9, fix_beta = TRUE)
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
colors <- rainbow(40)
clusters <- c()
for (i in c(1:40))
  clusters <- c(clusters, rep(i, 10))
results <- components(net)$membership
acc <- sum(results == clusters) / nrow(Y)
print(acc)
V(net)$cluster <- clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("face-recognition-sgl.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
