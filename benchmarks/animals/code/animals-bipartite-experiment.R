library(spectralGraphTopology)
library(igraph)
library(pals)
library(latex2exp)
library(extrafont)
library(viridis)
set.seed(0)

df <- read.csv("animals.txt", header = FALSE)
names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
n <- ncol(Y)
graph <- learn_bipartite_graph(t(Y), w0 = "qp", nu = 1e5, maxiter = 5e4, z = 1, abstol = 0)
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
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

xlab <- TeX("$\\mathit{p}$")
gr = .5 * (1 + sqrt(5))
eigvals <- c(eigenvalues(graph$Adjacency))
setEPS()
postscript("eigenvalues.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:length(eigvals)), eigvals, ylim=c(min(eigvals), max(eigvals)), xlab = xlab, ylab = "Eigenvalues of the Adjacency matrix",
     type = "p", col = "black", bg = inferno(3)[2], pch = 23)
grid()
dev.off()
embed_fonts("eigenvalues.ps", outfile="eigenvalues.ps")
