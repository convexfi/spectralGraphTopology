---
layout: default
title: Animals
parent: Clustering
nav_order: 1
---

# Animals
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Animal dataset

In this animal dataset (B. M. Lake an J. B. Tenenbaum 2010), there are 102 features which are
binary answers to questions such as "Does it has feathers?", "Is it warm-blooded?",
etc.  There are in total 33 animals to be clustered.


``` r
library(viridis)
library(spectralGraphTopology)
library(igraph)

df <- read.csv("animals.txt", header = FALSE)
names <- matrix(unlist(read.csv("animals_names.txt", header = FALSE)))
Y <- matrix(as.numeric(unlist(df)), nrow = nrow(df))
n <- nrow(Y)
graph <- learn_k_component_graph(cov(t(Y)) + diag(1/3, n, n), w0 = "qp",
                                 beta = 1, k = 10, verbose = FALSE)
net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
colors <- viridis(50, begin = 0, end = 1, direction = -1)
c_scale <- colorRamp(colors)
E(net)$color = apply(c_scale(abs(E(net)$weight) / max(abs(E(net)$weight))), 1,
                     function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
V(net)$color = "pink"
plot(net, vertex.label = names,
     vertex.size = 4,
     vertex.label.dist = 1,
     vertex.label.family = "Helvetica",
     vertex.label.cex = .8,
     vertex.label.color = "black")
```

![](animals_files/figure-markdown_github/unnamed-chunk-1-1.png)

- Lake, Brendan and Joshua Tenenbaum. "Discovering Structure by Learning Sparse Graphs."
  Proceedings of the 32nd Annual Meeting of the Cognitive Science Society CogSci 2010,
  Portland, Oregon, United States, 11-14 August, 2010, Cognitive Science Society, Inc., 2010. pp. 778-784.
