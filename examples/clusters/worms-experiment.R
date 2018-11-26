library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)
# number of points per cluster
N <- 100
worms <- shapes.worms(N)
graph <- learnGraphTopology(crossprod(t(worms$data)) + diag(rep(1/3, 2 * N)),
                            w0 = "naive", K = 2, beta = .25, ftol = 5e-4)
Lw <- graph$Lw
Lw <- diag(diag(Lw)) - Lw
net <- graph_from_adjacency_matrix(Lw, mode = "undirected", weighted = TRUE)
colors <- c("#706FD3", "#FF5252", "#33D9B2")
V(net)$cluster <- worms$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#00000'))
V(net)$color <- c(colors[1], colors[2])[worms$clusters]
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("worms.ps", family = "Times", height = 5, width = gr * 3.5)
#plot(worms$data, col=c(colors[5], colors[3])[worms$clusters], pch = 19,
#     xlab = "X", ylab = "Y")
plot(net, layout = worms$data, vertex.label = NA, vertex.size = 3)
dev.off()
colors <- c("#706FD3", "#FF5252", "#33D9B2")
setEPS()
postscript("worms_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(graph$elapsed_time, graph$loglike, type = "b", lty = 1, pch = 15, cex=.75, col = colors[1],
     xlab = "CPU time [seconds]", ylab = "")
grid()
lines(graph$elapsed_time, graph$obj_fun, type = "b", xaxt = "n", lty = 2, pch=16, cex=.75, col = colors[2])
lines(graph$elapsed_time, graph$obj_fun - graph$loglike, type = "b", xaxt = "n", lty = 3, pch=17, cex=.75,
      col = colors[3])
legend("topright", legend = c("likelihood", "posterior", "prior"),
       col=colors, pch=c(15, 16, 17), lty=c(1, 2, 3), bty="n")
dev.off()
embed_fonts("worms_trend.ps", outfile="worms_trend.ps")
