library(spectralGraphTopology)
library(igraph)
library(pals)
library(kernlab)
library(corrplot)
library(huge)
library(R.matlab)
data(iris)

n <- nrow(iris)
Y <- t(matrix(as.numeric(unlist(iris[, 1:4])), nrow = nrow(iris)))
names <- c(as.character(iris[, 5]))
unique_names <- c(as.character(unique(iris[, 5])))
S <- cov(Y)
## GGL
print("Connecting to MATLAB...")
matlab <- Matlab(port = 9999)
open(matlab)
print("success!")
A_mask <- matrix(1, nrow(iris), nrow(iris)) - diag(nrow(iris))
setVariable(matlab, A_mask = A_mask)
alpha <- 0.5
setVariable(matlab, alpha = alpha)
setVariable(matlab, S = S)
evaluate(matlab, "[Lggl,~,~] = estimate_ggl(S, A_mask, alpha, 1e-4, 1e-6, 1000, 1)")
Lggl <- getVariable(matlab, "Lggl")
Lw <- Lggl$Lggl
W <- diag(diag(Lw)) - Lw
ggl_net <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
###
colors <- c("#34495E", "#706FD3", "#FF5252")
clusters <- array(0, length(names))
for (i in c(1:length(unique_names)))
  clusters[unique_names[i] == names] <- i
V(ggl_net)$cluster <- clusters
E(ggl_net)$color <- apply(as.data.frame(get.edgelist(ggl_net)), 1,
                       function(x) ifelse(V(ggl_net)$cluster[x[1]] == V(ggl_net)$cluster[x[2]],
                                          colors[V(ggl_net)$cluster[x[1]]], brewer.greys(5)[2]))
V(ggl_net)$color <- colors[clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("iris_ggl.ps", family = "Times", height = 5, width = gr * 3.5)
plot(ggl_net, vertex.label = NA, vertex.size = 3)
dev.off()
