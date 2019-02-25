library(spectralGraphTopology)
library(igraph)
library(pals)
library(kernlab)
library(corrplot)
data(iris)

n <- nrow(iris)
Y <- t(matrix(as.numeric(unlist(iris[, 1:4])), nrow = nrow(iris)))
names <- c(as.character(iris[, 5]))
unique_names <- c(as.character(unique(iris[, 5])))
## Spectral Clustering
spec <- specc(t(Y), centers = 3)
spec_cluster <- as.vector(as.numeric(factor(spec)))
A <- matrix(0, n, n)
for (i in c(1:(n-1)))
  for (j in c((i+1):n))
    if (spec_cluster[i] == spec_cluster[j])
      A[i, j] <- 1
A <- A + t(A)
spec_net <- graph_from_adjacency_matrix(0.01 * A, mode = "undirected", weighted = TRUE)
colors <- c("#34495E", "#706FD3", "#FF5252")
### True clusters
clusters <- array(0, length(names))
for (i in c(1:length(unique_names)))
  clusters[unique_names[i] == names] <- i
###
V(spec_net)$cluster <- spec_cluster
V(spec_net)$color <- colors[clusters]
E(spec_net)$color <- "white"
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("iris_spec.ps", family = "Times", height = 5, width = gr * 3.5)
plot(spec_net, vertex.label = NA, vertex.size = 3)
dev.off()

