library(spectralGraphTopology)
library(igraph)
library(extrafont)

set.seed(1)

generate_helix <- function(n_points, n_turns, z_origin) {
  t <- seq(from = 0, to = n_turns * pi, by = n_turns * pi / n_points)[1:n_points]
  x <- cos(t)
  y <- sin(t)
  z <- z_origin + t
  return(list(x = x, y = y, z = z))
}

n_points <- 100
n_nodes <- 200
h1 <- generate_helix(n_points, 4, 0)
h2 <- generate_helix(n_points, 4, 1)
h1_sample_x <- rnorm(n_nodes, mean = h1$x, sd = 5e-2)
h1_sample_y <- rnorm(n_nodes, mean = h1$y, sd = 5e-2)
h1_sample_z <- rnorm(n_nodes, mean = h1$z, sd = 5e-2)
h2_sample_x <- rnorm(n_nodes, mean = h2$x, sd = 5e-2)
h2_sample_y <- rnorm(n_nodes, mean = h2$y, sd = 5e-2)
h2_sample_z <- rnorm(n_nodes, mean = h2$z, sd = 5e-2)

# learn underlying graph
Y <- rbind(cbind(h1_sample_x, h1_sample_y, h1_sample_z),
           cbind(h2_sample_x, h2_sample_y, h2_sample_z))
graph <- learn_k_component_graph(crossprod(t(Y)), k = 2, beta = 1, fix_beta = FALSE, abstol = 1e-3)
# construct network
net <- graph_from_adjacency_matrix(graph$adjacency, mode = "undirected", weighted = TRUE)
# use pretty colors for the nodes and edges
colors <- c("#706FD3", "#FF5252")
clusters <- c(rep(1, n_nodes), rep(2, n_nodes))
V(net)$cluster <- clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                     function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#000000'))
V(net)$color <- colors[clusters]
# set up network plot
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/helix3d.ps", height = 5, width = gr * 3.5)
plot(net, layout = Y[,2:3], vertex.label = NA, vertex.size = 3)
dev.off()
