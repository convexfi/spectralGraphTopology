[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spectralGraphTopology)](https://cran.r-project.org/package=spectralGraphTopology)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/spectralGraphTopology)](https://cran.r-project.org/package=spectralGraphTopology)
![CRAN Downloads Total](https://cranlogs.r-pkg.org/badges/grand-total/spectralGraphTopology?color=brightgreen)

[![codecov](https://codecov.io/gh/mirca/spectralGraphTopology/branch/master/graph/badge.svg)](https://codecov.io/gh/mirca/spectralGraphTopology)
[![Travis-CI-Badge](https://travis-ci.org/mirca/spectralGraphTopology.svg?branch=master)](https://travis-ci.org/mirca/spectralGraphTopology)
[![Build status](https://ci.appveyor.com/api/projects/status/vr62ddvc9xoabnwy?svg=true)](https://ci.appveyor.com/project/mirca/spectralgraphtopology-j05c9)
[![CircleCI](https://circleci.com/gh/mirca/spectralGraphTopology.svg?style=svg)](https://circleci.com/gh/mirca/spectralGraphTopology)
[![Docker Build Status](https://img.shields.io/docker/cloud/build/mirca/spectralgraphtopology.svg)](https://hub.docker.com/r/mirca/spectralgraphtopology/)
[![Build Status](https://dev.azure.com/jvmirca/spectralGraphTopology/_apis/build/status/mirca.spectralGraphTopology?branchName=master)](https://dev.azure.com/jvmirca/spectralGraphTopology/_build/latest?definitionId=1&branchName=master)

<p align="center">
  <img width = "200" src="./man/figures//circles3_reduced.gif"/>
</p>

**spectralGraphTopology** provides estimators to learn k-component, bipartite,
and k-component bipartite graphs from data by imposing spectral constraints
on the eigenvalues and eigenvectors of the Laplacian and adjacency matrices.
Those estimators leverages spectral properties of the graphical models as a
prior information, which turn out to play key roles in unsupervised machine
learning tasks such as community detection.

## Installation

From inside an R session, type:

```r
> devtools::install_github("dppalomar/spectralGraphTopology")
```

Alternatively, you can install the development version from GitHub:
```
$ git clone https://github.com/dppalomar/spectralGraphTopology.git
$ cd spectralGraphTopology
$ make build && make install
```

#### Microsoft Windows
On MS Windows environments, make sure to install the most recent version of ``Rtools``.

#### macOS
**spectralGraphTopology** depends on [`RcppArmadillo`](https://github.com/RcppCore/RcppArmadillo) which requires [`gfortran`](https://cloud.r-project.org/bin/macosx/tools/).

# Usage: clustering
We illustrate the usage of the package with simulated data, as follows:

```r
library(spectralGraphTopology)
library(clusterSim)
library(igraph)
set.seed(42)

# generate graph and data
n <- 50  # number of nodes per cluster
twomoon <- clusterSim::shapes.two.moon(n)  # generate datapoints
k <- 2  # number of components

# estimate underlying graph
S <- crossprod(t(twomoon$data))
graph <- learn_k_component_graph(S, k = k, beta = .5, verbose = FALSE, abstol = 1e-3)

# plot
# build network
net <- igraph::graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
# colorify nodes and edges
colors <- c("#706FD3", "#FF5252")
V(net)$cluster <- twomoon$clusters
E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
                      function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
                                        colors[V(net)$cluster[x[1]]], '#000000'))
V(net)$color <- colors[twomoon$clusters]
# plot nodes
plot(net, layout = twomoon$data, vertex.label = NA, vertex.size = 3)
```

<img src="man/figures/README-plot_k_component-1.png" width="75%" style="display: block; margin: auto;" />

For more examples, check out our [gallery](https://mirca.github.io/spectralGraphTopology).

## Contributing
We welcome all sorts of contributions. Please feel free to open an issue
to report a bug or discuss a feature request.

## Citation
If you made use of this software please consider citing:
* S. Kumar, J. Ying, J. V. de Miranda Cardoso, and D. P. Palomar (2019). A unified framework
  for structured graph learning via spectral constraints. https://arxiv.org/abs/1904.09792

In case you made use of the function `cluster_k_component_graph`, consider citing:
* N., Feiping, W., Xiaoqian, J., Michael I., and H., Heng. (2016).
  The Constrained Laplacian Rank Algorithm for Graph-based Clustering,
  AAAI'16. http://dl.acm.org/citation.cfm?id=3016100.3016174

## Links
Package: [GitHub](https://github.com/dppalomar/spectralGraphTopology)

README file: [GitHub-readme](https://raw.githack.com/dppalomar/spectralGraphTopology/master/README.html)
