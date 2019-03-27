library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

eps <- 1e-2
n_realizations <- 500
ratios <- c(10)
n <- 16
k <- 4
P <- diag(1, k)
mgraph <- sample_sbm(n, pref.matrix = P, block.sizes = c(rep(n / k, k)))
maxiter <- 5e4
relative_error_list <- list()
fscore_list <- list()
nll_list <- list()
objfun_list <- list()
constr_list <- list()

for (j in 1:length(ratios)) {
  t <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", t, "samples per node, t/n = ", ratios[j], "\n")
  for (r in 1:n_realizations) {
    print(r)
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(t, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_laplacian_matrix(S, w0 = "naive", k = 4, beta = 1e2, fix_beta = TRUE,
                                    edge_tol = eps, maxiter = maxiter,
                                    record_weights = TRUE, record_objective = TRUE)
    niter <- length(graph$loglike)
    relative_error <- array(0, niter)
    fscore <- array(0, niter)
    for (i in 1:niter) {
      Lw <- L(as.array(graph$w_seq[[i]]))
      relative_error[i] <- relativeError(Ltrue, Lw)
      fscore[i] <- metrics(Ltrue, Lw, eps)[1]
    }
    relative_error_list <- rlist::list.append(relative_error_list, relative_error)
    fscore_list <- rlist::list.append(fscore_list, fscore)
    nll_list <- rlist::list.append(nll_list, graph$loglike)
    objfun_list <- rlist::list.append(objfun_list, graph$obj_fun)
    constr_list <- rlist::list.append(constr_list, graph$obj_fun - graph$loglike)
  }
}

saveRDS(relative_error_list, file = "relerr_beta1e2.RDS")
saveRDS(fscore_list, file = "fscore_beta1e2.RDS")
saveRDS(nll_list, file = "nll_beta1e2.RDS")
saveRDS(objfun_list, file = "objfun_beta1e2.RDS")
saveRDS(constr_list, file = "constr_beta1e2.RDS")
