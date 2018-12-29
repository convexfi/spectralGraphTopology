library(spectralGraphTopology)
library(igraph)
library(pals)
library(R.matlab)
set.seed(0)

Nnodes <- 70

print("Connecting to MATLAB...")
matlab <- Matlab(port = 9998)
open(matlab)
print("success!")
A_mask <- matrix(1, Nnodes, Nnodes) - diag(Nnodes)
setVariable(matlab, A_mask = A_mask)

df <- read.csv("data.csv", header = FALSE, nrows = Nnodes)
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
Y <- Y[2:nrow(Y), 1:Nnodes]
df_names <- read.csv("labels.csv", header = FALSE, nrows = Nnodes)
names <- t(matrix(unlist(df_names), nrow = nrow(df_names)))
names <- names[2, 1:Nnodes]

S <- cov(Y)
setVariable(matlab, alpha = 1e-1)
setVariable(matlab, S = S)
evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-4, 1e-6, 200, 1)")
Lcgl <- getVariable(matlab, "Lcgl")
Lw <- Lcgl$Lcgl
W <- diag(diag(Lw)) - Lw

net <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
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
setEPS()
gr = .5 * (1 + sqrt(5))
postscript("../latex/figures/cancer-rna-graph-subset-ggl.ps", family = "Helvetica", height = 5, width = gr * 3.5)
#layout <- layout_in_circle(net, order = V(net))
#plot(net, layout = layout, vertex.label = names, vertex.size = 3)
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .4,
     vertex.label.color = "black")
dev.off()
