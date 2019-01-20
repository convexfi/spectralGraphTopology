library(spectralGraphTopology)
library(igraph)
library(pals)

Nnodes <- 70
df <- read.csv("data.csv", header = FALSE, nrows = Nnodes)
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
Y <- Y[2:nrow(Y), 1:Nnodes]
df_names <- read.csv("labels.csv", header = FALSE, nrows = Nnodes)
names <- t(matrix(unlist(df_names), nrow = nrow(df_names)))
names <- names[2, 1:Nnodes]

graph <- learnLaplacianGraphTopology(cov(Y), K = 5, w0 = "qp", beta = 5, Lwtol = 1e-6, maxiter = 100000)
print(graph$convergence)
net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
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
V(net)$color <- colors[clusters]
setEPS()
gr = .5 * (1 + sqrt(5))
postscript("../latex/figures/cancer-rna-graph-subset.ps", family = "Helvetica", height = 5, width = gr * 3.5)
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .4,
     vertex.label.color = "black")
dev.off()
