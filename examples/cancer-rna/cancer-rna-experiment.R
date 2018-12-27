library(spectralGraphTopology)
library(igraph)
library(pals)
set.seed(0)

Nnodes <- 80
df <- read.csv("data.csv", header = FALSE)
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
Y <- Y[2:nrow(Y), 1:Nnodes]
print(ncol(Y))
print(nrow(Y))
df_names <- read.csv("labels.csv", header = FALSE)
names <- t(matrix(unlist(df_names), nrow = nrow(df_names)))
names <- names[2, 1:Nnodes]
print(names)

N <- ncol(Y)
graph <- learnGraphTopology(cov(Y), K = 5, w0 = "naive", beta = 5, maxiter = 100000)
print(graph$lambda)
print(graph$convergence)
net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
##colors <- viridis(5, begin = 0, end = .1, direction = -1)
colors <- brewer.greys(10)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
setEPS()
postscript("cancer-rna-graph.ps")
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
