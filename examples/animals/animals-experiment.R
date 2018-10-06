library(spectralGraphTopology)
library(igraph)
library(pals)
library(viridis)
library(latex2exp)
set.seed(0)

df <- read.csv("animals.txt", header = FALSE)
names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
N <- ncol(Y)
graph <- learnGraphTopology(cov(Y) + diag(rep(1/3, N)), K = 10, beta = .5)

net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
#colors <- viridis(5, begin = 0, end = .1, direction = -1)
colors <- brewer.greys(3)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
#layout <- layout_in_circle(net, order = V(net))
#plot(net, layout = layout, vertex.label = NA, vertex.size = 3)
setEPS()
postscript("animals_graph_k10_alpha0.ps")
plot(net, vertex.label = names, vertex.size = 3,
     vertex.label.dist = 1, vertex.degree = pi/2,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black",
     main = TeX("$K = 10,\\alpha = 0$"))
dev.off()
#plot(graph$elapsed_time, graph$obj_fun, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
#      xlab = "CPU time [seconds]", ylab = "Objective function")
#plot(graph$elapsed_time, graph$loglike, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
#     xlab = "CPU time [seconds]", ylab = "Negative Log Likelihood")
