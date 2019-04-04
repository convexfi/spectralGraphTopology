library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
set.seed(42)

eps <- 1e-2
n_realizations <- 50
ratios <- c(30)
n <- 64
n1 <- 40
n2 <- 24
bipartite <- sample_bipartite(n1, n2, type="Gnp", p = .35, directed=FALSE)
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
    E(bipartite)$weight <- runif(gsize(bipartite), min = 1e-1, max = 3)
    #E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = .4)
    Ltrue <- as.matrix(laplacian_matrix(bipartite))
    #Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    #Y <- MASS::mvrnorm(t, mu = rep(0, n), Sigma = MASS::ginv(Ltrue + Lerdo))
    Y <- MASS::mvrnorm(t, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_bipartite_graph(S, z = abs(n2 - n1), w0 = "naive",
                                   nu = 1e1, abstol = 0, record_weights = TRUE)
    print(graph$convergence)
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
