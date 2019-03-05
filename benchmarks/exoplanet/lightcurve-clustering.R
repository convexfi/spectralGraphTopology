library(spectralGraphTopology)
library(corrplot)
library(pals)
library(igraph)
library(huge)
set.seed(42)

days <- 90 # days
obs_time <- 30 # minutes / per_observation
n_obs <- round(days / (obs_time / (60 * 24))) # unit of number of observations
obs <- c(1:n_obs)
days_per_observation <- days / n_obs
hours_per_observation <- days * 24 / n_obs

transits <- rep(0, n_obs)
noise_level = 0.01 ^ 2
n_curves <- 128
Y <- matrix(0, n_curves, n_obs)
# generate light curves with planets
p <- .4
meta_data <- matrix(0, n_curves, 6)
for (i in c(1:n_curves)) {
  noise <- MASS::mvrnorm(n_obs, mu = 0, Sigma = noise_level)
  if (runif(1) < p) {
    period <- round(runif(1, min = 5, max = 55) / days_per_observation)
    phase <- round(runif(1, min = 5, max = 35) / days_per_observation)
    duration <- round(runif(1, min = 8, max = 24) / hours_per_observation)
    depth <- runif(1, min = 0.001, max = 0.3)
    transits[(abs(obs - phase + period) %% period) < duration] <- depth
    meta_data[i, ] <- c(i, period, phase, duration, depth, 1)
    Y[i, ] <- 1 - transits + noise
  }
  else
    Y[i, ] <- 1 + noise
}

print(meta_data)
graph <- learn_laplacian_matrix(Y, k = 2, w0 = "qp", beta = 1e2)
print(graph$convergence)
graph_clr <- constr_laplacian_rank(Y, k = 2)
gr = .5 * (1 + sqrt(5))
est_net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
lasso_net <- graph_from_adjacency_matrix(graph_clr$Adjacency, mode = "undirected", weighted = TRUE)
colors <- brewer.blues(10)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(est_net)$color = c("red", "green")[meta_data[, 6] + 1]
setEPS()
postscript("est_graph.ps", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
E(lasso_net)$color = apply(c_scale(abs(E(lasso_net)$weight / max(E(lasso_net)$weight))), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(lasso_net)$color = c("red", "green")[meta_data[, 6] + 1]
setEPS()
postscript("lasso_graph.ps", height = 5, width = gr * 3.5)
plot(lasso_net, vertex.label = NA, vertex.size = 3)
dev.off()

