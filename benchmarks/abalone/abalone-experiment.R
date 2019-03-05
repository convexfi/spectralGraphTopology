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
unique_rings <- c(unique(rings))
clean_df <- abalone_df_ext[, !(names(abalone_df_ext) %in% drops)]
data_df <- data.matrix(clean_df)
# estimate graph
Y <- t(data_df)
k <- length(unique(rings))
graph <- constr_laplacian_rank(t(Y), k = k)
#graph <- learn_laplacian_matrix(S / max(S), w0 = "qp", k = k, beta = 5)
# plots
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
clusters <- array(0, length(rings))
for (i in c(1:length(unique_rings)))
  clusters[unique_rings[i] == rings] <- i
V(net)$cluster <- clusters
colors <- rainbow(length(unique_rings))
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("abalone-clr.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
