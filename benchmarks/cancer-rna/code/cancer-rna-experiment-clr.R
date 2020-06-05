library(spectralGraphTopology)
library(igraph)
library(pals)
set.seed(0)

Nnodes <- 801
df <- read.csv("data.csv", header = FALSE, nrows = Nnodes)
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
Y <- Y[2:nrow(Y), 1:Nnodes]
df_names <- read.csv("labels.csv", header = FALSE, nrows = Nnodes)
names <- t(matrix(unlist(df_names), nrow = nrow(df_names)))
names <- names[2, 1:Nnodes]

N <- ncol(Y)
graph <- cluster_k_component_graph(t(Y), k = 5, m = 7)
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
colors <- c("#34495E", "#706FD3", "#FF5252", "#33D9B2", "#34ACE0")
clusters <- array(0, length(names))
for (i in c(1:length(names))) {
  if (names[i] == "BRCA") {
    clusters[i] = 1
  } else if (names[i] == "COAD") {
    clusters[i] = 2
  } else if (names[i] == "LUAD") {
    clusters[i] = 3
  } else if (names[i] == "PRAD") {
    clusters[i] = 4
  } else {
    clusters[i] = 5
  }
}
V(net)$cluster <- clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], brewer.greys(5)[2]))
V(net)$color <- c(colors[1], colors[2], colors[3], colors[4], colors[5])[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("cancer-rna-graph-full-clr.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, vertex.label = NA,
     vertex.size = 3)
dev.off()
