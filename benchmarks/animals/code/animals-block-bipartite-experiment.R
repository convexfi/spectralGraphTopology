library(spectralGraphTopology)
library(igraph)
library(pals)
library(viridis)
library(latex2exp)
library(corrplot)
set.seed(0)

df <- read.csv("animals.txt", header = FALSE)
names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
n <- ncol(Y)
graph <- learn_bipartite_k_component_graph(t(Y), k = 5, w0 = "qp", nu = 1e5, beta = 1e6, fix_beta = TRUE,
                                       maxiter = 5e4, z = 3, edge_tol = 0, abstol = 0)
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
colors <- brewer.reds(100)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
setEPS()
postscript("../latex/figures/animals_graph_block_bipartite.ps")
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black")
dev.off()
