library(spectralGraphTopology)
library(corrplot)
library(pals)
library(extrafont)
library(igraph)
library(R.matlab)
set.seed(42)

eps <- 1e-1
n1 <- 40
n2 <- 24
n <- n1 + n2
pc <- .6

n_realizations <- 10
ratios <- c(10, 100, 500, 1000, 5000)

rel_err_sgl <- matrix(0, n_realizations, length(ratios))
rel_err_cgl <- matrix(0, n_realizations, length(ratios))
rel_err_naive <- matrix(0, n_realizations, length(ratios))
rel_err_qp <- matrix(0, n_realizations, length(ratios))
rel_err_pavez <- matrix(0, n_realizations, length(ratios))
fscore_sgl <- matrix(0, n_realizations, length(ratios))
fscore_cgl <- matrix(0, n_realizations, length(ratios))
fscore_naive <- matrix(0, n_realizations, length(ratios))
fscore_qp <- matrix(0, n_realizations, length(ratios))
fscore_pavez <- matrix(0, n_realizations, length(ratios))
recall_sgl <- matrix(0, n_realizations, length(ratios))
recall_cgl <- matrix(0, n_realizations, length(ratios))
recall_naive <- matrix(0, n_realizations, length(ratios))
recall_qp <- matrix(0, n_realizations, length(ratios))
recall_pavez <- matrix(0, n_realizations, length(ratios))
specificity_sgl <- matrix(0, n_realizations, length(ratios))
specificity_cgl <- matrix(0, n_realizations, length(ratios))
specificity_naive <- matrix(0, n_realizations, length(ratios))
specificity_qp <- matrix(0, n_realizations, length(ratios))
specificity_pavez <- matrix(0, n_realizations, length(ratios))
accuracy_sgl <- matrix(0, n_realizations, length(ratios))
accuracy_cgl <- matrix(0, n_realizations, length(ratios))
accuracy_naive <- matrix(0, n_realizations, length(ratios))
accuracy_qp <- matrix(0, n_realizations, length(ratios))
accuracy_pavez <- matrix(0, n_realizations, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab(port=9998)
is_matlab_open <- open(matlab)
cat("MATLAB connection status: ", is_matlab_open)
A_mask <- matrix(1, n, n) - diag(n)
setVariable(matlab, A_mask = A_mask)

n_ratios <- c(1:length(ratios))
for (j in n_ratios) {
  p <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", p, "samples per node, p/n = ", ratios[j], "\n")
  for (r in 1:n_realizations) {
    bipartite <- sample_bipartite(n1, n2, type="Gnp", p = pc, directed=FALSE)
    # randomly assign edge weights to connected nodes
    E(bipartite)$weight <- runif(gsize(bipartite), min = 1, max = 3)
    # get true Laplacian and Adjacency
    Ltrue <- as.matrix(laplacian_matrix(bipartite))
    Atrue <- diag(diag(Ltrue)) - Ltrue
    # set number of samples
    Y <- MASS::mvrnorm(p, rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    setVariable(matlab, S = S)
    setVariable(matlab, rr = diag(1/sqrt(diag(S))))
    evaluate(matlab, "[~, ~, Abi] = best_bipartite_approx(rr * S * rr)")
    Abi <- as.matrix(getVariable(matlab, "Abi")$Abi)
    setVariable(matlab, Abi = Abi)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(n)/p), 0)
    alpha_tmp <- alphas[1]
    rel_cgl <- Inf
    rel_pavez <- Inf
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-4, 1e-4, 100, 1)")
      evaluate(matlab, "[LcglA,~,~] = estimate_cgl(S, Abi, alpha, 1e-4, 1e-4, 100, 1)")
      Lcgl <- getVariable(matlab, "Lcgl")$Lcgl
      LcglA <- getVariable(matlab, "LcglA")$LcglA
      if (anyNA(Lcgl) || anyNA(LcglA)) {
        next
      }
      Acgl <- diag(diag(Lcgl)) - Lcgl
      AcglA <- diag(diag(LcglA)) - LcglA
      tmp_rel_cgl <- relativeError(Atrue, Acgl)
      tmp_rel_pavez <- relativeError(Atrue, AcglA)
      if (tmp_rel_cgl < rel_cgl) {
        alpha_tmp <- alpha
        rel_cgl <- tmp_rel_cgl
        metrics_cgl <- metrics(Atrue, Acgl, eps)
      }
      if (tmp_rel_pavez < rel_pavez) {
        rel_pavez <- tmp_rel_pavez
        metrics_pavez <- metrics(Atrue, AcglA, eps)
      }
    }
    Sinv <- MASS::ginv(S)
    w_naive <- spectralGraphTopology:::w_init(w0 = "naive", Sinv)
    w_qp <- spectralGraphTopology:::w_init(w0 = "qp", Sinv)
    graph <- learn_bipartite_graph(S, z = abs(n1 - n2), w0 = w_qp,
                                   nu = 1e5, abstol = 1e-4, maxiter = 1e5)
    print(graph$convergence)
    Anaive <- A(w_naive)
    Aqp <- A(w_qp)

    metrics_sgl <- metrics(Atrue, graph$Adjacency, eps)
    metrics_naive <- metrics(Atrue, Anaive, eps)
    metrics_qp <- metrics(Atrue, Aqp, eps)
    print(metrics_sgl)
    print(metrics_qp)
    print(metrics_pavez)
    print(metrics_cgl)

    rel_err_sgl[r, j] <- .5 * relativeError(Atrue, graph$Adjacency)
    rel_err_cgl[r, j] <- .5 * rel_cgl
    rel_err_naive[r, j] <- .5 * relativeError(Atrue, Anaive)
    rel_err_qp[r, j] <- .5 * relativeError(Atrue, Aqp)
    rel_err_pavez[r, j] <- .5 * rel_pavez
    print(rel_err_sgl[r, j])
    print(rel_err_qp[r, j])
    print(rel_pavez)
    print(rel_cgl)

    fscore_sgl[r, j] <-   metrics_sgl[1]
    fscore_naive[r, j] <- metrics_naive[1]
    fscore_qp[r, j] <-    metrics_qp[1]
    fscore_cgl[r, j] <-   metrics_cgl[1]
    fscore_pavez[r, j] <- metrics_pavez[1]

    recall_sgl[r, j] <-   metrics_sgl[2]
    recall_naive[r, j] <- metrics_naive[2]
    recall_qp[r, j] <-    metrics_qp[2]
    recall_cgl[r, j] <-   metrics_cgl[2]
    recall_pavez[r, j] <- metrics_pavez[2]

    specificity_sgl[r, j] <-   metrics_sgl[3]
    specificity_naive[r, j] <- metrics_naive[3]
    specificity_qp[r, j] <-    metrics_qp[3]
    specificity_cgl[r, j] <-   metrics_cgl[3]
    specificity_pavez[r, j] <- metrics_pavez[3]

    accuracy_sgl[r, j] <-   metrics_sgl[4]
    accuracy_naive[r, j] <- metrics_naive[4]
    accuracy_qp[r, j] <-    metrics_qp[4]
    accuracy_cgl[r, j] <-   metrics_cgl[4]
    accuracy_pavez[r, j] <- metrics_pavez[4]
  }
}

saveRDS(rel_err_sgl, file = "rel_err_sgl.rds")
saveRDS(rel_err_cgl, file = "rel_err_cgl.rds")
saveRDS(rel_err_naive, file = "rel_err_naive.rds")
saveRDS(rel_err_qp, file = "rel_err_qp.rds")
saveRDS(rel_err_pavez, file = "rel_err_pavez.rds")

saveRDS(fscore_sgl, file = "fscore_sgl.rds")
saveRDS(fscore_naive, file = "fscore_naive.rds")
saveRDS(fscore_qp, file = "fscore_qp.rds")
saveRDS(fscore_cgl, file = "fscore_cgl.rds")
saveRDS(fscore_pavez, file = "fscore_pavez.rds")

saveRDS(recall_sgl, file = "recal_sgl.rds")
saveRDS(recall_naive, file = "recal_naive.rds")
saveRDS(recall_qp, file = "recal_qp.rds")
saveRDS(recall_cgl, file = "recal_cgl.rds")
saveRDS(recall_pavez, file = "recal_pavez.rds")

saveRDS(specificity_sgl, file = "specificity_sgl.rds")
saveRDS(specificity_naive, file = "specificity_naive.rds")
saveRDS(specificity_qp, file = "specificity_qp.rds")
saveRDS(specificity_cgl, file = "specificity_cgl.rds")
saveRDS(specificity_pavez, file = "specificity_pavez.rds")

saveRDS(accuracy_sgl, file = "accuracy_sgl.rds")
saveRDS(accuracy_naive, file = "accuracy_naive.rds")
saveRDS(accuracy_qp, file = "accuracy_qp.rds")
saveRDS(accuracy_cgl, file = "accuracy_cgl.rds")
saveRDS(accuracy_pavez, file = "accuracy_pavez.rds")

