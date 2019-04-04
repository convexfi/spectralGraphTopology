library(spectralGraphTopology)
library(igraph)
library(pals)

synthetic_df <- read.csv("synthetic_control.data", header = FALSE, sep = " ")
Y <- matrix(as.numeric(unlist(synthetic_df)), nrow = nrow(synthetic_df))
Y <- scale(Y)
print(dim(Y))
labels <- c()
unique_labels <- c(1:6)
print(length(unique_labels))
for (i in unique_labels)
  labels <- c(labels, rep(i, 100))
graph <- constr_laplacian_rank(Y, k = length(unique_labels))
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected",
                                   weighted = TRUE)
V(net)$cluster <- labels
colors <- c("#574b90", "#f78fb3", "#3dc1d3",
            "#546de5", "#e15f41", "#f5cd79")
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[labels]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("synthetic-control.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
