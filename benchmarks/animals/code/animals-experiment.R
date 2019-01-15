library(spectralGraphTopology)
library(igraph)
library(pals)
library(latex2exp)
library(huge)
set.seed(0)

run_animals <- function(K) {
  df <- read.csv("animals.txt", header = FALSE)
  names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
  Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
  N <- ncol(Y)
  graph <- learnLaplacianGraphTopology(cov(Y) + diag(rep(1/3, N)), w0 = "qp",
                                       K = K, beta = .5, alpha = 1e-1)
  net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
  #colors <- viridis(5, begin = 0, end = .1, direction = -1)
  colors <- brewer.greys(10)
  c_scale <- colorRamp(colors)
  E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                       function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
  V(net)$color = "pink"
  setEPS()
  postscript(paste0("../latex/figures/animals_graph_k", toString(K), ".ps"))
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

for(k in c(1, 2, 4, 8, 10)) {
  print(paste("running animals exp. for K =", toString(k)))
  run_animals(K = k)
}
