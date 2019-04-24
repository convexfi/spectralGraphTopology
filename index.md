---
layout: default
title: Home
nav_order: 1
description: "*spectralGraphTopology* is a software tool to estimate "
permalink: /
---

# Learning similarities in graph data
{: .fs-9 }

**spectralGraphTopology** provides estimators of the Laplacian and adjacency matrices of
graphs by leveraing prior information in their structural form.

{: .fs-6 .fw-300 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/dppalomar/spectralGraphTopology){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Getting started

### Dependencies

The R version of **spectralGraphTopology** is build on top of awesome R packages including **Rcpp**,
**RcppEigen**, and **osqp**. All these packages can be installed via CRAN.

### Installation

2. The **development** version can be installed via GitHub as
```bash
devtools::install_github("dppalomar/spectralGraphTopology")
```
<small>You must have previously installed the **devtools** package.</small>


### Tutorials

- [See the package vignette](http://mirca.github.io/spectralGraphTopology/_static/getting_started.html) for a
detailed description of the mathematical methods that are available in **spectralGraphTopology**.

---

## About the project

**spectralGraphTopology** is developed on [GitHub](http://github.com/dppalomar/spectralGraphTopology)
by [Ze Vinicius](http://mirca.github.io), [Daniel Palomar](http://www.danielppalomar.com), Jiaxi Ying
and Sandeep Kumar.

### License

**spectralGraphTopology** is distributed by an
[GPL 3.0 License](https://github.com/dppalomar/spectralGraphTopology/blob/master/LICENSE).

### Contributing

We welcome all sorts of contributions. Please feel free to open an issue to report a bug or discuss a feature request in [our GitHub repo](https://github.com/dppalomar/spectralGraphTopology).

### Citation

If this package has been useful to you in any way, give us a star on [GitHub](http://github.com/dppalomar/spectralGraphTopology) :)
Additionally, if you've used **spectralGraphTopology** on your research, please consider citing the following resource:

- S. Kumar, J. Ying, J. V. de M. Cardoso and D. P. Palomar (2019). A unified framework for structured graph learning
  via spectral constraints. https://arxiv.org/pdf/1904.09792.pdf
