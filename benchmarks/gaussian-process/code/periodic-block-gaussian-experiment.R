library(igraph)
library(corrplot)
library(spectralGraphTopology)
library(pals)
library(huge)
library(extrafont)
library(R.matlab)
set.seed(0)

n <- 128
p <- 30 * n
S1 <- periodic_kernel(.5, sqrt(2), 10, n/4) + gaussian_kernel(.5, 5e-1, n/4)
S2 <- periodic_kernel(.5, sqrt(2), 10, n/4) + gaussian_kernel(.5, 5e-1, n/4)
S3 <- periodic_kernel(.5, sqrt(2), 10, n/4) + gaussian_kernel(.5, 5e-1, n/4)
S4 <- periodic_kernel(.5, sqrt(2), 10, n/4) + gaussian_kernel(.5, 5e-1, n/4)
Strue <- blockDiag(S1, S2, S3, S4) + wiener_kernel(n)
Strue <- Strue / mean(Strue)
Ptrue <- spectralGraphTopology:::inv_sympd(Strue)
Wtrue <- diag(diag(Ptrue)) - Ptrue
true_net <- graph_from_adjacency_matrix(abs(Wtrue), mode = "undirected", weighted = TRUE)
Y <- MASS::mvrnorm(p, mu = rep(0, n), Sigma = Strue)
SCM <- cov(Y)
# SGL
graph <- learn_laplacian_matrix(SCM, k = 1, w0 = "qp")
# GLasso
glasso <- huge(SCM, method = "glasso")
Pglasso <- matrix(glasso$icov[[10]], nrow = 128)
# GGL
print("Connecting to MATLAB...")
matlab <- Matlab()
open(matlab)
print("success!")
A_mask <- matrix(1, n, n) - diag(n)
setVariable(matlab, A_mask = A_mask)
setVariable(matlab, SCM = SCM)
alpha = 1e-2
setVariable(matlab, alpha = alpha)
evaluate(matlab, "[Lggl,~,~] = estimate_ggl(SCM, A_mask, alpha, 1e-6, 1e-6, 40, 1)")
Lggl <- getVariable(matlab, "Lggl")$Lggl
# construct nets
glasso_net <- graph_from_adjacency_matrix(abs(diag(diag(Ptrue)) - Ptrue), mode = "undirected", weighted = TRUE)
est_net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
blues <- brewer.blues(100)
reds <- brewer.reds(100)
blue_scale <- colorRamp(blues)
red_scale <- colorRamp(reds)
E(est_net)$color = apply(blue_scale(E(est_net)$weight / max(E(est_net)$weight)),
                         1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
non_zero_edges <- Ainv(Wtrue)[Ainv(Wtrue) != 0]
pos_edges <- non_zero_edges > 0
neg_edges <- !pos_edges
edge_colors <- (blue_scale(pos_edges * E(true_net)$weight / max(pos_edges * E(true_net)$weight)) +
                red_scale(neg_edges * abs(E(true_net)$weight) / max(neg_edges * abs(E(true_net)$weight))))
edge_colors <- 255 * edge_colors / max(edge_colors)
E(true_net)$color = apply(edge_colors, 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
print(graph$convergence)

gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/periodic_gaussian_block_sgl_graph.ps",
           family = "Times", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_block_true_graph.ps",
           family = "Times", height = 5, width = gr * 3.5)
plot(true_net, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_block_laplacian.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Laplacian / max(graph$Laplacian), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_block_ggl.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(Lggl / max(Lggl), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_block_covariance.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(Strue / max(Strue), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_block_precision.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(Ptrue / max(Ptrue), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
setEPS()
postscript("../latex/figures/periodic_gaussian_block_glasso_precision.ps",
           family = "Times", height = 5, width = gr * 3.5)
corrplot(Pglasso / max(Pglasso), is.corr = FALSE, method = "square",
         addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()
#colors <- c("#0B032D", "#843B62", "#F67E7D", "#6ABA81")
#pch <- c(15, 7, 8, 9)
#setEPS()
#postscript("../latex/figures/sample_periodic_gaussian_block.ps",
#           family = "Times")#, height = 5, width = gr * 7)
#plot(c(1:n), Y[1,], type = "l", lty = 1, col = colors[1],
#     xlab = "Time stamp", ylab = "Signal strength")
##lines(c(1:p), Y[,2], type = "l", lty = 1, col = colors[2])
##lines(c(1:p), Y[,3], type = "l", lty = 1, col = colors[3])
##lines(c(1:p), Y[,4], type = "l", lty = 1, col = colors[4])
#dev.off()

