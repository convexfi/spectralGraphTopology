library(fastDummies)
library(spectralGraphTopology)
library(igraph)
library(pals)

yeast_df <- read.csv("yeast.data", sep = " ", header = FALSE)
Y <- matrix(as.numeric(unlist(yeast_df[, 2:8])), nrow = nrow(yeast_df))
names <- c(as.character(yeast_df[, 9]))
unique_names <- c(as.character(unique(names)))
#graph <- learn_laplacian_matrix(Y, m = 30, w0 = "qp", k = 3, beta = 1e2)
#print(graph$beta_seq)
#print(graph$convergence)
graph <- constr_laplacian_rank(Y, k = length(unique_names))
saveRDS(graph$Adjacency, "clr-adjacency.rds")
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
clusters <- array(0, length(names))
for (i in c(1:length(unique_names)))
  clusters[unique_names[i] == names] <- i
V(net)$cluster <- clusters
colors <- rainbow(length(unique_names))
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("yeast-clr.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA, vertex.size = 3)
dev.off()
