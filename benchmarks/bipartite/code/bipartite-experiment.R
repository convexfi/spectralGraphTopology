library(spectralGraphTopology)
library(corrplot)
library(pals)
library(extrafont)
library(igraph)
library(R.matlab)
set.seed(234)

n1 <- 40
n2 <- 24
n <- n1 + n2
pc <- .6

n_realizations <- 20
ratios <- c(10, 30, 100, 250, 500, 1000, 2000, 4000, 8000)
rel_err_sgl <- matrix(0, n_realizations, length(ratios))
rel_err_cgl <- matrix(0, n_realizations, length(ratios))
rel_err_cglA <- matrix(0, n_realizations, length(ratios))
rel_err_naive <- matrix(0, n_realizations, length(ratios))
rel_err_qp <- matrix(0, n_realizations, length(ratios))
fscore_sgl <- matrix(0, n_realizations, length(ratios))
fscore_cgl <- matrix(0, n_realizations, length(ratios))
fscore_cglA <- matrix(0, n_realizations, length(ratios))
fscore_naive <- matrix(0, n_realizations, length(ratios))
fscore_qp <- matrix(0, n_realizations, length(ratios))

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
    w_true <- Linv(Ltrue)
    # set number of samples
    Y <- MASS::mvrnorm(p, rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    setVariable(matlab, S = S)
    evaluate(matlab, "[~, ~, A] = best_bipartite_approx(S)")
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(n)/p), 0)
    rel_cgl <- Inf
    rel_cglA <- Inf
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-6, 1e-6, 300, 1)")
      Lcgl <- getVariable(matlab, "Lcgl")
      Acgl <- diag(diag(Lcgl$Lcgl)) - Lcgl$Lcgl
      if (anyNA(Lcgl$Lcgl)) {
        next
      }
      tmp_rel_cgl <- relativeError(Atrue, Acgl)
      if (tmp_rel_cgl < rel_cgl) {
        rel_cgl <- tmp_rel_cgl
        fs_cgl <- Fscore(Atrue, Acgl, 1e-1)
      }
      evaluate(matlab, "[LcglA,~,~] = estimate_cgl(S, A, alpha, 1e-6, 1e-6, 300, 1)")
      LcglA <- getVariable(matlab, "LcglA")
      AcglA <- diag(diag(LcglA$LcglA)) - LcglA$LcglA
      if (anyNA(LcglA$LcglA)) {
        next
      }
      tmp_rel_cglA <- relativeError(Atrue, AcglA)
      if (tmp_rel_cglA < rel_cglA) {
        rel_cglA <- tmp_rel_cglA
        fs_cglA <- Fscore(Atrue, AcglA, 1e-1)
      }
    }
    Sinv <- MASS::ginv(S)
    w_naive <- spectralGraphTopology:::w_init(w0 = "naive", Sinv)
    w_qp <- spectralGraphTopology:::w_init(w0 = "qp", Sinv)
    graph <- learn_bipartite_graph(S, z = abs(n2 - n1), w0 = w_qp, maxiter = 1e5)
    print(graph$convergence)
    Anaive <- A(w_naive)
    Aqp <- A(w_qp)
    rel_sgl = relativeError(Atrue, graph$Adjacency)
    fs_sgl = Fscore(Atrue, graph$Adjacency, 1e-1)
    rel_naive = relativeError(Atrue, Anaive)
    fs_naive = Fscore(Atrue, Anaive, 1e-1)
    rel_qp = relativeError(Atrue, Aqp)
    fs_qp = Fscore(Atrue, Aqp, 1e-1)
    rel_err_sgl[r, j] <- rel_sgl
    fscore_sgl[r, j] <- fs_sgl
    rel_err_cgl[r, j] <- rel_cgl
    fscore_cgl[r, j] <- fs_cgl
    rel_err_cglA[r, j] <- rel_cglA
    fscore_cglA[r, j] <- fs_cglA
    rel_err_naive[r, j] <- rel_naive
    fscore_naive[r, j] <- fs_naive
    rel_err_qp[r, j] <- rel_qp
    fscore_qp[r, j] <- fs_qp
  }
}
saveRDS(rel_err_sgl, file = "rel-err-SGL.rds")
saveRDS(fscore_sgl, file = "fscore-SGL.rds")
saveRDS(rel_err_cgl, file = "rel-err-CGL.rds")
saveRDS(fscore_cgl, file = "fscore-CGL.rds")
saveRDS(rel_err_cglA, file = "rel-err-CGLA.rds")
saveRDS(fscore_cglA, file = "fscore-CGLA.rds")
saveRDS(rel_err_naive, file = "rel-err-naive.rds")
saveRDS(fscore_naive, file = "fscore-naive.rds")
saveRDS(rel_err_qp, file = "rel-err-QP.rds")
saveRDS(fscore_qp, file = "fscore-QP.rds")

## build the network
#net <- graph_from_adjacency_matrix(Aw, mode = "undirected", weighted = TRUE)
#V(net)$type = V(bipartite)$type
## plot network
#colors <- brewer.blues(20)
#c_scale <- colorRamp(colors)
#E(net)$color = apply(c_scale(E(net)$weight / max(E(net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
#E(bipartite)$color = apply(c_scale(E(bipartite)$weight / max(E(bipartite)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
#gr = .5 * (1 + sqrt(5))
##setEPS()
##postscript("test.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
##plot(net, vertex.label = NA,
##     vertex.color=c("#706FD3", "#33D9B2")[V(net)$type+1],
##     vertex.size = 3)
##dev.off()
##embed_fonts("test.ps", outfile="test.ps")
##
##setEPS()
##postscript("test_true.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
##plot(bipartite, vertex.label = NA,
##     vertex.color=c("#706FD3", "#33D9B2")[V(bipartite)$type+1],
##     vertex.size = 3)
##dev.off()
##embed_fonts("test_true.ps", outfile="test_true.ps")
#
#setEPS()
#postscript("est_mat.ps", family = "Times", height = 5, width = gr * 3.5)
#corrplot(graph$Aw / max(graph$Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
#dev.off()
#setEPS()
#postscript("true_mat.ps", family = "Times", height = 5, width = gr * 3.5)
#corrplot(Aw / max(Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
#dev.off()
#
## plot convergence trend
#niter <- length(graph$loglike)
#iterations <- c(1:niter)
#setEPS()
#postscript("test_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(iterations, graph$obj_fun, type = "b", lty = 1, pch = 15, cex=.75, col = "#706FD3",
#     xlab = "Iteration number", ylab = "Objective value")
#grid()
#legend("topright", legend = c("objective-function"),
#       col=c("#706FD3"), pch=c(15), lty=c(1), bty="n")
#dev.off()
#embed_fonts("test_trend.ps", outfile="test_trend.ps")
