---
layout: default
title: Bipartite
parent: Structured graphs
nav_order: 1
---

# Bipartite
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

``` r
library(spectralGraphTopology)
library(igraph)
library(viridis)
library(corrplot)
n1 <- 10
n2 <- 6
n <- n1 + n2
pc <- .9

bipartite <- sample_bipartite(n1, n2, type="Gnp", p = pc, directed=FALSE)
# randomly assign edge weights to connected nodes
E(bipartite)$weight <- runif(gsize(bipartite), min = 0, max = 1)
# Erdo-Renyi as noise model
a <- .35
erdos_renyi <- erdos.renyi.game(n, p = .35)
E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
# get true Laplacian and Adjacency
Ltrue <- as.matrix(laplacian_matrix(bipartite))
Atrue <- diag(diag(Ltrue)) - Ltrue
Lnoisy = Ltrue + Lerdo
Anoisy <- diag(diag(Lnoisy)) - Lnoisy
# set number of samples
Y <- MASS::mvrnorm(10 * n, rep(0, n), Sigma = MASS::ginv(Lnoisy))
# compute sample covariance matrix
S <- cov(Y)
# estimate Adjacency matrix
graph <- learn_bipartite_graph(S, w0 = "qp", verbose = FALSE)
# Plot Adjacency matrices: true, noisy, and estimated
corrplot(Atrue / max(Atrue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
corrplot(Anoisy / max(Anoisy), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
corrplot(graph$Adjacency / max(graph$Adjacency), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
# Plot networks
estimated_bipartite <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
V(estimated_bipartite)$type <- c(rep(0, 10), rep(1, 6))
la = layout_as_bipartite(estimated_bipartite)
colors <- viridis(20, begin = 0, end = 1, direction = -1)
c_scale <- colorRamp(colors)
E(estimated_bipartite)$color = apply(c_scale(E(estimated_bipartite)$weight / max(E(estimated_bipartite)$weight)), 1,
                          function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
E(bipartite)$color = apply(c_scale(E(bipartite)$weight / max(E(bipartite)$weight)), 1,
                      function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
la = la[, c(2, 1)]
# Plot networks: true and estimated
plot(bipartite, layout = la, vertex.color=c("red","black")[V(bipartite)$type + 1],
     vertex.shape = c("square", "circle")[V(bipartite)$type + 1],
     vertex.label = NA, vertex.size = 3)
plot(estimated_bipartite, layout = la, vertex.color=c("red","black")[V(estimated_bipartite)$type + 1],
     vertex.shape = c("square", "circle")[V(estimated_bipartite)$type + 1],
     vertex.label = NA, vertex.size = 3)
```

True adjacency             |  Noisy adjacency          |  Estimated adjacency
:-------------------------:|:-------------------------:|:------------------------:
![](bipartite_files/figure-markdown_github/unnamed-chunk-1-1.png) | ![](bipartite_files/figure-markdown_github/unnamed-chunk-1-2.png) | ![](bipartite_files/figure-markdown_github/unnamed-chunk-1-3.png)

True graph                 |  Estimated graph
:-------------------------:|:-------------------------:
![](bipartite_files/figure-markdown_github/unnamed-chunk-1-4.png) | ![](bipartite_files/figure-markdown_github/unnamed-chunk-1-5.png)

