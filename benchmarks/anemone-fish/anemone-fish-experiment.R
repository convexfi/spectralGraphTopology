library(spectralGraphTopology)
library(igraph)
library(corrplot)

fish_names <-c("clarkii" ,"akindynos" ,"chrysopterus" ,"chrysogaster" ,"bicinctus",
               "tricinctus" ,"perideraion" ,"allardi" ,"melanopus" ,"ocellaris" ,"percula",
               "leucokranos","polymnus","omanensis","rubrocinctus","sandaracinos","akallopisos",
               "ephippium","frenatus","fuscocaudatus","latezonatus","latifasciatus","mccullochi",
               "nigripes", "sebae","biaculeatus")
A <- as.matrix(read.csv("anemonefish.txt"))
k <- nrow(A)
l <- ncol(A)
S <- rbind(cbind(matrix(0, k, k), A), cbind(t(A), matrix(0, l, l)))

graph <- learn_adjacency_and_laplacian(S, k = 5, w0 = "qp", beta = 1e3, fix_beta = TRUE, alpha = 5e-1, z = 18, maxiter = 1e5, edge_tol = 0)
corrplot(graph$Adjacency / max(graph$Adjacency), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
