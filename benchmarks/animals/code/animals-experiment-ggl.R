library(spectralGraphTopology)
library(igraph)
library(pals)
library(latex2exp)

set.seed(42)

df <- read.csv("ggl-laplacian-alpha005.txt", header = FALSE)
names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
laplacian <- matrix(as.numeric(unlist(df)), nrow = nrow(df))
adj <- diag(diag(laplacian)) - laplacian

net <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
colors <- brewer.reds(100)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
#layout <- layout_in_circle(net, order = V(net))
#plot(net, layout = layout, vertex.label = NA, vertex.size = 3)
setEPS()
postscript("../latex/figures/animals_graph_ggl_alpha01.ps")
plot(net, vertex.label = names, vertex.size = 3,
     vertex.label.dist = 1, vertex.degree = pi/2,
     vertex.label.family = "Helvectica",
     vertex.label.cex = .8,
     vertex.label.color = "black")
dev.off()
#plot(graph$elapsed_time, graph$obj_fun, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
#      xlab = "CPU time [seconds]", ylab = "Objective function")
#plot(graph$elapsed_time, graph$loglike, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
#     xlab = "CPU time [seconds]", ylab = "Negative Log Likelihood")
