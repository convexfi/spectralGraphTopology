library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(R.matlab)
set.seed(123)

eps <- 1e-2
n_realizations <- 100
ratios <- c(100)
n <- 64
grid <- make_lattice(length = sqrt(n), dim = 2)
relative_error_list <- list()
fscore_list <- list()
nll_list <- list()
objfun_list <- list()
constr_list <- list()
time_list <- list()

print("Connecting to MATLAB...")
matlab <- Matlab(port=9999)
open(matlab)
print("success!")
A_mask <- matrix(1, 64, 64) - diag(64)
setVariable(matlab, A_mask = A_mask)

for (j in 1:length(ratios)) {
  t <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", t, "samples per node, t/n = ", ratios[j], "\n")
  for (r in 1:n_realizations) {
    E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(grid))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(t, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    setVariable(matlab, S = S)
    evaluate(matlab, "[~,Laps,~,convergence] = estimate_cgl_list(S, A_mask, 0, 1e-4, 1e-4, 100, 1)")
    Laplacian_matrices <- getVariable(matlab, "Laps")$Laps
    times <- getVariable(matlab, "convergence")$convergence
    times <- cumsum(times)
    niter <- length(times)
    relative_error <- array(0, niter)
    fscore <- array(0, niter)
    nll <- array(0, niter)
    for (i in 1:niter) {
      Lw <- matrix(unlist(Laplacian_matrices[[i]]), n, n)
      relative_error[i] <- relativeError(Ltrue, Lw)
      fscore[i] <- metrics(Ltrue, Lw, eps)[1]
      eigvals <- c(eigenvalues(Lw))
      nll[i] <- spectralGraphTopology:::laplacian.likelihood(Lw, eigvals[2:length(eigvals)], S)
    }
    print(relative_error[1])
    print(relative_error[niter])
    print(fscore[1])
    print(fscore[niter])
    relative_error_list <- rlist::list.append(relative_error_list, relative_error)
    fscore_list <- rlist::list.append(fscore_list, fscore)
    nll_list <- rlist::list.append(nll_list, nll)
    #objfun_list <- rlist::list.append(objfun_list, graph$obj_fun)
    #constr_list <- rlist::list.append(constr_list, graph$obj_fun - graph$loglike)
    time_list <- rlist::list.append(time_list, times)
  }
}

saveRDS(relative_error_list, file = "relerr-cgl.RDS")
saveRDS(fscore_list, file = "fscore-cgl.RDS")
saveRDS(nll_list, file = "nll-cgl.RDS")
#saveRDS(objfun_list, file = "objfun-cgl.RDS")
#saveRDS(constr_list, file = "constr-cgl.RDS")
saveRDS(time_list, file = "time-cgl.RDS")
