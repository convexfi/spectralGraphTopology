library(spectralGraphTopology)
library(igraph)
library(pals)
library(matrixStats)
library(viridis)
library(latex2exp)
set.seed(0)

nnodes <- 1000
# get feature data
df <- read.csv("Hymapfeat.txt", header = FALSE, sep = " ")
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
ns <- sample(c(1:nrow(df)), nnodes)
Y <- Y[, ns]
Y <- scale(Y)
# get labels
df_names <- read.csv("Hymaplabel.txt", header = FALSE)
names <- t(matrix(unlist(df_names), nrow = nrow(df_names)))
names <- names[, ns]
k <- length(unique(names))
print(k)
# estimate graph
S <- cov(Y)
graph <- learn_laplacian_matrix(S/max(S), w0 = "naive", k = k, beta = 1.5, alpha = 1.5, maxiter = 1e5)
print(graph$convergence)
print(graph$lambda)
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
# colorify edges and nodes
colors <- c("#B33771", "#1B9CFC", "#EAB543", "#58B19F", "#D6A2E8", "#BDC581", "#182C61")
V(net)$cluster <- names + 1
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                      function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                         colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[names + 1]
setEPS()
postscript("hyperspectral_graph.ps")
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
