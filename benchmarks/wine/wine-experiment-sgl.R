library(igraph)
library(spectralGraphTopology)
library(pals)

wine_df <- read.csv("wine.data", header = FALSE)
labels <- array(unlist(wine_df[, 1]))
unique_labels <- c(unique(labels))
# estimate graph
Y <- matrix(unlist(wine_df[, 2:ncol(wine_df)]), nrow = nrow(wine_df))
k <- length(unique(labels))
Y <- scale(Y)
graph <- learn_k_component_graph(Y, k = k, beta = 2, fix_beta = TRUE, maxiter = 1e5)
# plots
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
clusters <- array(0, length(labels))
for (i in c(1:k))
  clusters[unique_labels[i] == labels] <- i
results <- components(net)$membership
acc <- sum(results == clusters) / length(labels)
print(acc)
V(net)$cluster <- clusters
colors <- c("#fab1a0", "#a29bfe", "#55efc4")
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("wine-sgl.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
