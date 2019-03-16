library(igraph)
library(viridis)
library(spectralGraphTopology)
library(pals)

# get data
k2c2_df <- read.csv("K2Campaign2_lightcurves_NoMotionDetrending.csv")
n_features <- ncol(k2c2_df)
ids <- k2c2_df[, n_features]
data <- data.matrix(k2c2_df[, 1:(n_features - 1)])
# get labels
results_df <- read.csv("K2Campaign2_lightcurves_NoMotionDetrending_classifications.csv")
labels <- results_df[, 2]
unique_labels <- unique(labels)
# estimate graph
k <- length(unique_labels)
S <- cov(t(data))
print(dim(S))
graph <- constr_laplacian_rank(data, k = k)
## plots
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
clusters <- array(0, length(labels))
for (i in c(1:length(unique_labels)))
  clusters[unique_labels[i] == labels] <- i
V(net)$cluster <- clusters
colors <- rainbow(length(unique_labels))
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("k2-clusters-clr.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
