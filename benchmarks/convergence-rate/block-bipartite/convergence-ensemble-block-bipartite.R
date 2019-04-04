library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
set.seed(42)

eps <- 5e-2
n_realizations <- 100
ratios <- c(100)

n <- 32
# randomly assign edge weights to connected nodes
b1 <- sample_bipartite(6, 4, type="Gnp", p = .8, directed=FALSE)
b2 <- sample_bipartite(7, 5, type="Gnp", p = .7, directed=FALSE)
b3 <- sample_bipartite(5, 5, type="Gnp", p = .9, directed=FALSE)

#erdos_renyi <- erdos.renyi.game(n, p = .35)
maxiter <- 5e4
relative_error_list <- list()
fscore_list <- list()
nll_list <- list()
objfun_list <- list()
constr_list <- list()
time_list <- list()

for (j in 1:length(ratios)) {
  t <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", t, "samples per node, t/n = ", ratios[j], "\n")
  for (r in 1:n_realizations) {
    E(b1)$weight <- runif(gsize(b1), min = .1, max = 3)
    E(b2)$weight <- runif(gsize(b2), min = .1, max = 3)
    E(b3)$weight <- runif(gsize(b3), min = .1, max = 3)
    Lw1 <- as.matrix(laplacian_matrix(b1))
    Lw2 <- as.matrix(laplacian_matrix(b2))
    Lw3 <- as.matrix(laplacian_matrix(b3))
    Ltrue <- blockDiag(Lw1, Lw2, Lw3)
    #Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    #Y <- MASS::mvrnorm(t, mu = rep(0, n), Sigma = MASS::ginv(Ltrue + Lerdo))
    Y <- MASS::mvrnorm(t, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_adjacency_and_laplacian(S, z = 0, w0 = "naive", k = 3, nu = 1e4, beta = 1e2,
                                           fix_beta = TRUE, maxiter = 1e3, reltol = 1e-3,
                                           record_weights = TRUE, record_objective = TRUE)
    niter <- length(graph$loglike)
    relative_error <- array(0, niter)
    fscore <- array(0, niter)
    for (i in 1:niter) {
      Lw <- L(as.array(graph$w_seq[[i]]))
      relative_error[i] <- relativeError(Ltrue, Lw)
      fscore[i] <- metrics(Ltrue, Lw, eps)[1]
    }
    print(relative_error[1])
    print(relative_error[niter])
    print(fscore[1])
    print(fscore[niter])
    relative_error_list <- rlist::list.append(relative_error_list, relative_error)
    fscore_list <- rlist::list.append(fscore_list, fscore)
    nll_list <- rlist::list.append(nll_list, graph$loglike)
    objfun_list <- rlist::list.append(objfun_list, graph$obj_fun)
    constr_list <- rlist::list.append(constr_list, graph$obj_fun - graph$loglike)
    time_list <- rlist::list.append(time_list, graph$elapsed_time)
  }
}

saveRDS(relative_error_list, file = "relerr.RDS")
saveRDS(fscore_list, file = "fscore.RDS")
saveRDS(nll_list, file = "nll.RDS")
saveRDS(objfun_list, file = "objfun.RDS")
saveRDS(constr_list, file = "constr.RDS")
saveRDS(time_list, file = "time.RDS")
