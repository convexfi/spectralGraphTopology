library(igraph)
library(pals)
library(latex2exp)
library(knor)
set.seed(0)

run_animals <- function(k) {
  df <- read.csv("animals.txt", header = FALSE)
  names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
  Y <- t(matrix(as.numeric(unlist(df)), nrow = nrow(df)))
  kmeans <- Kmeans(t(Y), centers = k)
  A <- matrix(0, 33, 33)
  for (i in c(1:32))
    for (j in c((i+1):33))
      if (kmeans$cluster[i] == kmeans$cluster[j])
        A[i, j] <- 1
  A <- A + t(A)
  net <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE)
  colors <- c("#EA2027", "#A3CB38", "#12CBC4", "#FDA7DF", "#D980FA",
              "#833471", "#9980FA", "#0652DD", "#1B1464", "#009432")
  V(net)$cluster <- kmeans$cluster
  V(net)$color <- colors[kmeans$cluster]
  E(net)$color <- "white"
  gr = .5 * (1 + sqrt(5))
  setEPS()
  postscript("../latex/figures/animals-kmeans.ps")
  plot(net, vertex.label = names,
       vertex.size = 4,
       vertex.label.dist = 1,
       vertex.label.family = "Helvetica",
       vertex.label.cex = .8,
       vertex.label.color = "black")
  dev.off()
}

for(k in c(10)) {
  print(paste("running animals exp. for k =", toString(k)))
  run_animals(k)
}
