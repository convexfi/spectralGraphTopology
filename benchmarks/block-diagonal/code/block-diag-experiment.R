library(igraph)
library(corrplot)
library(spectralGraphTopology)
library(pals)
library(extrafont)
library(huge)
library(R.matlab)
set.seed(0)

N <- 20
T <- 30 * N
K <- 4
P <- diag(1, K)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
Ltrue <- as.matrix(laplacian_matrix(mgraph))
Wtrue <- diag(diag(Ltrue)) - Ltrue
# Erdo-Renyi as noise model
p <- .35
a <- .45
erdos_renyi <- erdos.renyi.game(N, p)
E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
# Noisy Laplacian
Lnoisy <- Ltrue + Lerdo
Wnoisy <- diag(diag(Lnoisy)) - Lnoisy
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Lnoisy))
S <- cov(Y)
Sinv <- MASS::ginv(S)
w_qp <- spectralGraphTopology:::w_init("qp", Sinv)
w_naive <- spectralGraphTopology:::w_init("naive", Sinv)
Lqp <- L(w_qp)
Lnaive <- L(w_naive)
graph <- learn_laplacian_matrix(S, k = K, w0 = "naive", beta = 400, fix_beta = TRUE,
                                ub = 2*N, record_objective = TRUE, alpha = 0.1,
                                abstol = 0, maxiter = 1e4)
graph_glasso <- huge(S, method = "glasso", nlambda = 100)
Lglasso <- matrix(graph_glasso$icov[[100]], N, N)

matlab <- Matlab()
open(matlab)
A_mask <- matrix(1, N, N) - diag(N)
alpha = 5e-2
setVariable(matlab, S = S)
setVariable(matlab, A_mask = A_mask)
setVariable(matlab, alpha = alpha)
evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-4, 1e-4, 100, 1)")
Lcgl <- getVariable(matlab, "Lcgl")

est_net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
noisy_net <- graph_from_adjacency_matrix(Wnoisy, mode = "undirected", weighted = TRUE)
true_net <- graph_from_adjacency_matrix(Wtrue, mode = "undirected", weighted = TRUE)
colors <- brewer.blues(20)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(noisy_net)$color = apply(c_scale(E(noisy_net)$weight / max(E(noisy_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(true_net)$color = apply(c_scale(E(true_net)$weight / max(E(true_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))

cat("SGL results: \n")
print(relativeError(Ltrue, graph$Laplacian))
print(metrics(Ltrue, graph$Laplacian, 1e-3))
print(spectralGraphTopology:::laplacian.likelihood(graph$Laplacian, graph$lambda, S))
cat("QP results: \n")
print(relativeError(Ltrue, Lqp))
print(metrics(Ltrue, Lqp, 1e-3))
lambda_qp <- eigenvalues(Lqp)
print(spectralGraphTopology:::laplacian.likelihood(Lqp, lambda_qp[lambda_qp > 1e-9], S))
cat("Naive results: \n")
print(relativeError(Ltrue, Lnaive))
print(metrics(Ltrue, Lnaive, 1e-3))
lambda_naive <- eigenvalues(Lnaive)
print(spectralGraphTopology:::laplacian.likelihood(Lnaive, lambda_naive[lambda_naive > 1e-9], S))
cat("Glasso results: \n")
print(relativeError(Ltrue, Lglasso))
print(metrics(Ltrue, Lglasso, 1e-3))
lambda_glasso <- eigenvalues(Lglasso)
print(spectralGraphTopology:::laplacian.likelihood(Lglasso, lambda_glasso[lambda_glasso > 1e-9], S))
cat("CGL results: \n")
print(relativeError(Ltrue, Lcgl$Lcgl))
print(metrics(Ltrue, Lcgl$Lcgl, 1e-3))
lambda_cgl <- eigenvalues(Lcgl$Lcgl)
print(spectralGraphTopology:::laplacian.likelihood(Lcgl$Lcgl, lambda_cgl[lambda_cgl > 1e-9], S))


gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/est_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/noisy_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(noisy_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/true_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(true_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/est_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Adjacency / max(graph$Adjacency), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/noisy_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Wnoisy / max(Wnoisy), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/true_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(Wtrue / max(Wtrue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()

colors <- c("#706FD3", "#FF5252", "#33D9B2")
setEPS()
postscript("../latex/figures/bd_trend.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(c(1:length(graph$loglike)), graph$loglike, type = "b", lty = 1, pch = 15, cex=.75, col = colors[1],
     xlab = "Iteration Number", ylab = "", ylim=c(0, 15))
grid()
lines(c(1:length(graph$loglike)), graph$obj_fun, type = "b", xaxt = "n", lty = 2, pch=16, cex=.75, col = colors[2])
lines(c(1:length(graph$loglike)), graph$obj_fun - graph$loglike, type = "b", xaxt = "n", lty = 3, pch=17, cex=.75,
      col = colors[3])
legend("topright", legend = c("likelihood", "posterior", "prior"),
       col=colors, pch=c(15, 16, 17), lty=c(1, 2, 3), bty="n")
dev.off()
embed_fonts("../latex/figures/bd_trend.ps", outfile="../latex/figures/bd_trend.ps")
