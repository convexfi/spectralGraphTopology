library(spectralGraphTopology)
library(corrplot)
library(pals)
library(igraph)
set.seed(42)

days <- 90 # days
obs_time <- 30 # minutes / per_observation
n_obs <- round(days / (obs_time / (60 * 24))) # unit of number of observations
obs <- c(1:n_obs)
days_per_observation <- days / n_obs
hours_per_observation <- days * 24 / n_obs

transits <- rep(0, n_obs)
noise_level = 0.01 ^ 2
n_curves <- 70
Y <- matrix(0, n_curves, n_obs)
# generate light curves with planets
p <- .1
for (i in c(1:n_curves)) {
  period <- round(runif(1, min = 10, max = 30) / days_per_observation)
  phase <- round(runif(1, min = 5, max = 15) / days_per_observation)
  duration <- round(runif(1, min = 8, max = 16) / hours_per_observation)
  depth = runif(1, min = 0.001, 0.3)
  transits[(abs(obs - phase + period) %% period) < duration] <- depth
  noise <- MASS::mvrnorm(n_obs, mu = 0, Sigma = noise_level)
  if (runif(1) < p) {
    Y[i, ] <- 1 - transits + noise
    print(i)
  }
  else
    Y[i, ] <- 1 + noise
}

Y <- t(Y)
S <- cov(Y)
Sinv <- MASS::ginv(S)
graph <- learn_laplacian_matrix(S / max(S), k = 4, w0 = "qp", beta = .25, alpha = 1e-1)
print(graph$convergence)

gr = .5 * (1 + sqrt(5))
setEPS()
postscript("est_mat.ps", family = "Times", height = 5, width = gr * 3.5)
corrplot(graph$Aw / max(graph$Aw), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
dev.off()

est_net <- graph_from_adjacency_matrix(graph$Aw, mode = "undirected", weighted = TRUE)
colors <- brewer.blues(10)
c_scale <- colorRamp(colors)
E(est_net)$color = apply(c_scale(E(est_net)$weight / max(E(est_net)$weight)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(est_net)$color = "pink"
setEPS()
postscript("est_graph.ps", family = "Times", height = 5, width = gr * 3.5)
plot(est_net, vertex.label = NA, vertex.size = 3)
dev.off()
