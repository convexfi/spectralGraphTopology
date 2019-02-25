library(spectralGraphTopology)
library(igraph)
library(kernlab)
library(pals)
set.seed(0)

n <- 801
df <- read.csv("data.csv", header = FALSE, nrows = n)
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
Y <- Y[2:nrow(Y), 1:n]
df_names <- read.csv("labels.csv", header = FALSE, nrows = n)
names <- t(matrix(unlist(df_names), nrow = nrow(df_names)))
names <- names[2, 1:n]

spec <- specc(t(Y), centers = 5)
spec_cluster <- as.vector(as.numeric(factor(spec)))
A <- matrix(0, n, n)
for (i in c(1:(n-1)))
  for (j in c((i+1):n))
    if (spec_cluster[i] == spec_cluster[j])
      A[i, j] <- 1
A <- A + t(A)
spec_net <- graph_from_adjacency_matrix(0.05 * A, mode = "undirected", weighted = TRUE)
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
V(spec_net)$cluster <- spec_cluster
V(spec_net)$color <- colors[clusters]
E(spec_net)$color <- "white"
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/spec-cancer.ps", family = "Times", height = 5, width = gr * 3.5)
plot(spec_net, vertex.label = NA, vertex.size = 3)
dev.off()
