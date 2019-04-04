library(spectralGraphTopology)
library(igraph)
library(pals)
library(viridis)
library(latex2exp)
library(huge)
set.seed(0)

df <- read.csv("animals.txt", header = FALSE)
names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
N <- ncol(Y)
graph <- learn_dregular_graph(cov(Y) + diag(rep(1/3, N)), w0 = "qp",
                              k = 1, beta = 10, eta = 1e3, fix_beta = TRUE, maxiter = 1e5)
print(diag(graph$Laplacian))
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
colors <- brewer.reds(100)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
setEPS()
postscript(paste0("../latex/figures/animals_dregular_graph.ps"))
#layout <- layout_in_circle(net, order = V(net))
#plot(net, layout = layout, vertex.label = names, vertex.size = 3)
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black")
dev.off()
