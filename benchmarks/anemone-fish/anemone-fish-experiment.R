library(spectralGraphTopology)
library(igraph)

#fish_names <-c("clarkii" ,"akindynos" ,"chrysopterus" ,"chrysogaster" ,"bicinctus",
#               "tricinctus" ,"perideraion" ,"allardi" ,"melanopus" ,"ocellaris" ,"percula",
#               "leucokranos","polymnus","omanensis","rubrocinctus","sandaracinos","akallopisos",
#               "ephippium","frenatus","fuscocaudatus","latezonatus","latifasciatus","mccullochi",
#               "nigripes", "sebae","biaculeatus")
#A <- matrix(runif(6), nrow = 2, ncol = 3)
#k <- nrow(A)
#l <- ncol(A)
#Y <- rbind(cbind(matrix(0, k, k), A), cbind(t(A), matrix(0, l, l)))
#S <- cov(Y) + diag(1/3, 5)
#Sinv <- MASS::ginv(Y)
#print(Sinv)

n1 <- 10
n2 <- 4
n <- 3 * n1 + 3 * n2 - 10
p <- 10000 * n

# randomly assign edge weights to connected nodes
b1 <- sample_bipartite(n1, n2, type="Gnp", p = .7, directed=FALSE)
b2 <- sample_bipartite(n1-4, n2, type="Gnp", p = .8, directed=FALSE)
b3 <- sample_bipartite(n1-6, n2, type="Gnp", p = .9, directed=FALSE)
E(b1)$weight <- runif(gsize(b1), min = 1, max = 3)
E(b2)$weight <- runif(gsize(b2), min = 1, max = 3)
E(b3)$weight <- runif(gsize(b3), min = 1, max = 3)
Lw1 <- as.matrix(laplacian_matrix(b1))
Aw1 <- diag(diag(Lw1)) - Lw1
Lw2 <- as.matrix(laplacian_matrix(b2))
Aw2 <- diag(diag(Lw2)) - Lw2
Lw3 <- as.matrix(laplacian_matrix(b3))
Aw3 <- diag(diag(Lw3)) - Lw3
Lw <- blockDiag(Lw1, Lw2, Lw3)
Aw <- blockDiag(Aw1, Aw2, Aw3)

# set number of samples
Y <- MASS::mvrnorm(p, rep(0, n), Sigma = MASS::ginv(Lw))
S <- cov(Y)
Sinv <- MASS::ginv(S)
