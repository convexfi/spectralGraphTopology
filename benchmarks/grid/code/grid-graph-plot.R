library(igraph)
library(spectralGraphTopology)
library(corrplot)
library(scales)
library(viridis)
library(R.matlab)

set.seed(0)

N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
T <- as.integer(100 * N)
E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
Ltrue <- as.matrix(laplacian_matrix(grid))
Wtrue <- diag(diag(Ltrue)) - Ltrue
Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
S <- cov(Y)

print("Connecting to MATLAB...")
matlab <- Matlab(port = 9998)
open(matlab)
print("success!")
A_mask <- matrix(1, 64, 64) - diag(64)
setVariable(matlab, A_mask = A_mask)
setVariable(matlab, S = S)
alpha = 5e-3
setVariable(matlab, alpha = alpha)
evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-6, 1e-6, 40, 1)")
Lcgl <- getVariable(matlab, "Lcgl")
graph <- learn_k_component_graph(S, w0 = "qp", beta = 20, alpha = alpha, abstol = 1e-5, fix_beta = TRUE)

eps <- 5e-2
# compute adjacency matrices
graph$Adjacency[graph$Adjacency < eps] <- 0
W_cgl <- diag(diag(Lcgl$Lcgl)) - Lcgl$Lcgl
W_cgl[W_cgl < eps] <- 0
# compute grids
grid_spec <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
grid_cgl <- graph_from_adjacency_matrix(W_cgl, mode = "undirected", weighted = TRUE)

colors <- viridis(5, begin = 0, end = .75, direction = -1)
c_scale <- colorRamp(colors)
E(grid_spec)$color = apply(c_scale(E(grid_spec)$weight / max(E(grid_spec)$weight)), 1,
                          function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(grid_cgl)$color = apply(c_scale(E(grid_cgl)$weight / max(E(grid_cgl)$weight)), 1,
                          function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(grid)$color = apply(c_scale(E(grid)$weight / max(E(grid)$weight)), 1,
                      function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(grid_spec)$color = "grey"
V(grid_cgl)$color = "grey"
V(grid)$color = "grey"

la_spec <- layout_on_grid(grid_spec)
la_cgl <- layout_on_grid(grid_cgl)
la_true <- layout_on_grid(grid)

relativeError(Ltrue, graph$Laplacian)
metrics(Ltrue, graph$Laplacian, eps)
relativeError(Ltrue, Lcgl$Lcgl)
metrics(Ltrue, Lcgl$Lcgl, eps)

gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/true_grid.ps", height = 5, width = gr * 3.5)
plot(grid, layout = la_true, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/sgl_grid.ps", height = 5, width = gr * 3.5)
plot(grid_spec, layout = la_spec, vertex.label = NA, vertex.size = 3)
dev.off()
setEPS()
postscript("../latex/figures/cgl_grid.ps", height = 5, width = gr * 3.5)
plot(grid_cgl, layout = la_cgl, vertex.label = NA, vertex.size = 3)
dev.off()
