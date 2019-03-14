library(spectralGraphTopology)
library(igraph)
library(pals)
library(kernlab)
library(corrplot)
data(musk)

n_features <- ncol(musk) - 1
Y <- t(matrix(as.numeric(unlist(musk[, 1:n_features])), nrow = nrow(musk)))
labels <- c(as.numeric(musk[, n_features+1]))
unique_labels <- c(as.numeric(unique(musk[, n_features+1])))
graph <- learn_laplacian_matrix(t(Y), k = 2, w0 = "naive", beta = 1e2, m = 3, alpha = 1e-1, tol = 1e-6)
print(graph$convergence)
print(graph$beta_seq)
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
colors <- c("#34495E", "#706FD3")#, "#FF5252")#, "#33D9B2", "#34ACE0")
clusters <- array(0, length(labels))
for (i in c(1:length(unique_labels)))
  clusters[unique_labels[i] == labels] <- i
V(net)$cluster <- clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("musk.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
