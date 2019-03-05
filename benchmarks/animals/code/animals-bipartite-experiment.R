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
graph <- learn_adjacency_and_laplacian(cov(Y) + diag(rep(1/3, n)), w0 = "qp", k = 10, z = 1,
                                       beta1 = 5e-1, beta2 = 1e3, maxiter = 1e5)
print(graph$convergence)
net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
colors <- brewer.reds(100)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
setEPS()
postscript("../latex/figures/animals_graph_bipartite.ps")
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black")
dev.off()
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("est_animal_adjacency_matrix.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Aw / max(graph$Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()