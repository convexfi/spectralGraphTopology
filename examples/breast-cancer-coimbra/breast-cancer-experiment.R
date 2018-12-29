library(spectralGraphTopology)
library(igraph)
library(pals)
library(latex2exp)
library(huge)
library(matrixStats)
set.seed(0)

df <- read.csv("breast-cancer.data", header = FALSE)
#names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
names <- Y[nrow(Y),]
print(names)
Y <- Y[1:(nrow(Y)-1),]
mean_row <- rowMeans(Y)
std_row <- rowSds(Y)
Y <- (Y - mean_row) / std_row
N <- ncol(Y)
graph <- learnGraphTopology(crossprod(Y) + diag(1/4, N), K = 2, w0 = "naive", beta = 10, alpha = 1e-1)
print(graph$lambda)
print(graph$convergence)
graph$W[graph$W < 5e-2] = 0
net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
##colors <- viridis(5, begin = 0, end = .1, direction = -1)
colors <- brewer.greys(10)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
setEPS()
postscript("breast-cancer.ps")
#layout <- layout_in_circle(net, order = V(net))
#plot(net, layout = layout, vertex.label = names, vertex.size = 3)
plot(net, vertex.label = names,
     vertex.size = 3,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black")
dev.off()

n <- length(graph$loglike)
plot(c(1:n), graph$obj_fun - graph$loglike)
print(graph$obj_fun[n] - graph$loglike[n])
