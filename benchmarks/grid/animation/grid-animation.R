library(clusterSim)
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(extrafont)

set.seed(1)

N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
T <- as.integer(100 * N)
E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
Ltrue <- as.matrix(laplacian_matrix(grid))
Wtrue <- diag(diag(Ltrue)) - Ltrue
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
S <- cov(Y)

alpha = 1.5e-3
graph <- learn_k_component_graph(S, w0 = "qp", beta = 20, alpha = alpha, abstol = 1e-5, fix_beta = TRUE, record_weights = TRUE, record_objective = TRUE)

# pretty colors
colors <- c("#706FD3", "#FF5252", "#33D9B2")
gr = .5 * (1 + sqrt(5))
index <- c(1, 100, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, length(graph$obj_fun))
j <- 0
colors <- viridis(5, begin = 0, end = .75, direction = -1)
c_scale <- colorRamp(colors)
for(i in index) {
  print(j)
  net <- graph_from_adjacency_matrix(A(c(as.matrix(graph$w_seq[[i]]))), mode = "undirected", weighted = TRUE)
  layout <- layout_on_grid(net)
  V(net)$color = "grey"
  E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1,
                       function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
  filename <- paste0("test", j, ".pdf")
  pdf(filename, height = 5, width = 3.5 * gr)
  plot(net, layout = layout, vertex.label = NA, vertex.size = 4, edge.width = E(grid_spec)$weight)
  dev.off()
  j <- j + 1
}

system("convert -dispose Background -delay 40 *.pdf grid_reduced.gif")
file.remove(list.files(pattern=".pdf"))
