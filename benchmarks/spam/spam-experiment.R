library(spectralGraphTopology)
library(igraph)
library(pals)
library(kernlab)
library(corrplot)
data(spam)
set.seed(42)

n_nodes <- 128
nodes <- sample(1:dim(spam)[1], n_nodes)
n_features <- ncol(spam) - 1
Y <- t(matrix(as.numeric(unlist(spam[nodes, 1:n_features])), nrow = n_nodes))
labels <- c(as.character(spam[nodes, (n_features+1)]))
unique_labels <- c(as.character(unique(spam[nodes, (n_features+1)])))
print(unique_labels)
S <- crossprod(Y) / n_features
graph <- learn_k_component_graph(S, k = 2, w0 = "qp", beta = 1e5, maxiter = 1e5)
print(graph$convergence)
print(graph$lambda)
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
colors <- c("#34495E", "#FF5252")#, "#33D9B2", "#34ACE0")
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
postscript("spam.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
