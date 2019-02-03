library(fastDummies)
library(igraph)
library(viridis)
library(spectralGraphTopology)
library(pals)

abalone_df <- read.csv("abalone.data")
# use one-hot enconding for the "sex" categorical variable
abalone_df_ext <- dummy_cols(abalone_df)
# let's drop "sex" categorical variable and "ring", which is the variable
# we want to cluster for
drops <- c("sex", "rings")
rings <- matrix(unlist(abalone_df_ext["rings"]))
n <- 64
clean_df <- abalone_df_ext[, !(names(abalone_df_ext) %in% drops)]
data_df <- data.matrix(clean_df)
# estimate graph
Y <- t(data_df[1:n,])
S <- cov(Y)
k <- length(unique(rings[1:n]))
print(k)
graph <- learn_laplacian_matrix(S / max(S), w0 = "qp", k = k, beta = 5)
# plots
print(graph$convergence)
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
colors <- rainbow(length(unique(rings[1:n])))
V(net)$cluster <- rings[1:n]
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                      function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                         colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[c(1:length(rings[1:n]))]
plot(net, vertex.label = NA, vertex.size = 3)
