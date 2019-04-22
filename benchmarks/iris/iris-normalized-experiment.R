library(spectralGraphTopology)
library(igraph)
library(pals)
data(iris)

n <- nrow(iris)
Y <- matrix(as.numeric(unlist(iris[, 1:4])), nrow = nrow(iris))
names <- c(as.character(iris[, 5]))
unique_names <- c(as.character(unique(iris[, 5])))
S <- cov(t(Y))
graph <- learn_normalized_laplacian(S, k = 3, maxiter = 5000)
Adjacency <- diag(diag(graph$NormalizedLaplacian)) - graph$NormalizedLaplacian
net <- graph_from_adjacency_matrix(Adjacency, mode = "undirected", weighted = TRUE)
clusters <- array(0, length(names))
for (i in c(1:length(unique_names)))
  clusters[unique_names[i] == names] <- i
V(net)$cluster <- clusters
colors <- c("#34495E", "#706FD3", "#FF5252")
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("iris-normalized.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
