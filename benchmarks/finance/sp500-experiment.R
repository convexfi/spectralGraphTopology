library(huge)
library(pals)
library(spectralGraphTopology)
data(stockdata)

Y <- log(stockdata$data[2:1258,]/stockdata$data[1:1257,])
labels <- stockdata$info[, 2]
unique_labels <- unique(labels)
k <- length(unique_labels)
S <- cov(Y)
print(dim(S))
graph <- learn_laplacian_matrix(S/max(S), k = 1, beta = 10, maxiter = 1e5)
print(graph$convergence)
iter <- length(graph$obj_fun)
clusters <- numeric(k)
for (i in 1:k) {
  clusters[labels == unique_labels[i]] <- i
}
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
colors <- rainbow(k)
V(net)$cluster <- clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                      function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                         colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("finance.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
