library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)
# number of nodes per cluster
N <- 100
# generate datapoints
twomoon <- shapes.two.moon(N)
# estimate underlying graph
graph <- learnLaplacianGraphTopology(crossprod(t(twomoon$data)) + diag(rep(1/3, 2 * N)),
                            w0 = "naive", K = 2, beta = .25, ftol = 1e-3, Lwtol = 1e-3)
# build network
net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
# colorify nodes and edges
colors <- c("#706FD3", "#FF5252", "#33D9B2")
V(net)$cluster <- twomoon$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#00000'))
V(net)$color <- c(colors[1], colors[2])[twomoon$clusters]
# plot network
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/twomoon.ps", family = "Times", height = 5, width = gr * 3.5)
plot(net, layout = twomoon$data, vertex.label = NA, vertex.size = 3)
dev.off()
colors <- c("#706FD3", "#FF5252", "#33D9B2")
# plot convergence trend
setEPS()
postscript("../latex/figures/twomoon_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:length(graph$loglike)), graph$loglike, type = "b", lty = 1, pch = 15, cex=.75, col = colors[1],
     xlab = "Iteration Number", ylab = "")
grid()
lines(c(1:length(graph$loglike)), graph$obj_fun, type = "b", xaxt = "n", lty = 2, pch=16, cex=.75, col = colors[2])
lines(c(1:length(graph$loglike)), graph$obj_fun - graph$loglike, type = "b", xaxt = "n", lty = 3, pch=17, cex=.75,
      col = colors[3])
legend("topright", legend = c("likelihood", "posterior", "prior"),
       col=colors, pch=c(15, 16, 17), lty=c(1, 2, 3), bty="n")
dev.off()
embed_fonts("../latex/figures/twomoon_trend.ps", outfile="../latex/figures/twomoon_trend.ps")
