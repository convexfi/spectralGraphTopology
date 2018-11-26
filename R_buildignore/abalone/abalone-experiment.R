library(fastDummies)
library(igraph)
library(viridis)
library(spectralGraphTopology)

abalone_df <- read.csv("abalone.data")
# use one-hot enconding for the "sex" categorical variable
abalone_df_ext <- dummy_cols(abalone_df)
# let's drop "sex" categorical variable and "ring", which is the variable
# we want to cluster for
drops <- c("sex", "rings")
rings <- matrix(unlist(abalone_df_ext["rings"]))
clean_df <- abalone_df_ext[, !(names(abalone_df_ext) %in% drops)]
data_df <- data.matrix(clean_df)
# estimate graph
N <- 1000
graph <- learnGraphTopology(crossprod(t(data_df[1:N,])) + rep(1/3, N),
                            w0 = "naive", K = 29, beta = .25, ftol = 1e-3)
# plots
net <- graph_from_adjacency_matrix(graph$W, mode = "undirected", weighted = TRUE)
c_scale <- colorRamp(inferno(5, begin = 0, end = 1, direction = -1))
# Applying the color scale to edge weights.
# rgb method is to convert colors to a character vector.
E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = rings[1:N]
plot(net, vertex.label = NA, vertex.size = 3)#, edge.width=E(net)$weight)
plot(graph$elapsed_time, graph$obj_fun, type = "b", pch=19, cex=.6, col = scales::alpha("black", .5),
      xlab = "CPU time [seconds]", ylab = "Objective function")
