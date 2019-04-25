library(spectralGraphTopology)
library(igraph)
library(pals)
library(viridis)
library(latex2exp)
library(huge)
set.seed(0)

run_animals <- function(k) {
  df <- read.csv("animals.txt", header = FALSE)
  names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
  Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
  N <- ncol(Y)
  if (k == 1)
    graph <- learn_k_component_graph(cov(Y) + diag(1/3, N, N), w0 = "qp", k = k,
                                    beta = .5, alpha = 5e-2)
  else
    graph <- learn_laplacian_matrix(cov(Y) + diag(1/3, N, N), w0 = "qp", beta = 1, k = k)
  print(graph$elapsed_time)
  net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
  colors <- brewer.reds(100)
  c_scale <- colorRamp(colors)
  E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                       function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
  V(net)$color = "pink"
  setEPS()
  postscript(paste0("../latex/figures/animals_graph_k", toString(k), ".ps"))
  #layout <- layout_in_circle(net, order = V(net))
  #plot(net, layout = layout, vertex.label = names, vertex.size = 3)
  plot(net, vertex.label = names,
       vertex.size = 3,
       vertex.label.dist = 1,
       vertex.label.family = "Helvetica",
       vertex.label.cex = .8,
       vertex.label.color = "black")
  dev.off()
}

for(k in c(1, 10)) {
  print(paste("running animals exp. for K =", toString(k)))
  run_animals(k)
}
