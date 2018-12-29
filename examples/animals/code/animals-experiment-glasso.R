library(igraph)
library(pals)
library(latex2exp)
library(huge)
set.seed(0)

df <- read.csv("animals.txt", header = FALSE)
names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
N <- ncol(Y)
graph <- huge(cov(Y) + diag(rep(1/3, N)), method = "glasso")
print(graph$lambda)
W <- matrix(graph$icov[[10]], nrow = 33)
W <- diag(diag(W)) - W
net <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
colors <- brewer.greys(10)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
setEPS()
postscript(paste0("../latex/figures/animals_glasso.ps"))
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black")
dev.off()
