library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)

set.seed(123)
# number of points per cluster
N <- 50
twomoon <- shapes.two.moon(N)
graph <- learnGraphTopology(crossprod(t(twomoon$data)) + diag(rep(1/3, 2 * N)),
                            w0 = "naive", K = 2, beta = .25, ftol = 9e-4)
net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
colors <- inferno(5, begin = 0, end = 1, direction = -1)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = c(colors[5], colors[3])[twomoon$clusters]
plot(twomoon$data, col=c(colors[5], colors[3])[twomoon$clusters], pch = 19,
     xlab = "X", ylab = "Y")
plot(net, layout = twomoon$data, vertex.label = NA, vertex.size = 3)
plot(graph$elapsed_time, graph$obj_fun, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
     xlab = "CPU time [seconds]", ylab = "Objective function")
