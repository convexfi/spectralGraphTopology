#' Package spectralGraphTopology
#'
#' This package provides estimators to learn k-component, bipartite,
#' and k-component bipartite graphs from data by imposing spectral constraints
#' on the eigenvalues and eigenvectors of the Laplacian and adjacency matrices.
#' Those estimators leverages spectral properties of the graphical models as a
#' prior information, which turn out to play key roles in unsupervised machine
#' learning tasks such as community detection.
#'
#' @section Functions:
#' \code{\link{learn_k_component_graph}}
#' \code{\link{learn_bipartite_graph}}
#' \code{\link{learn_bipartite_k_component_graph}}
#' \code{\link{cluster_k_component_graph}}
#' \code{\link{L}}
#' \code{\link{A}}
#'
#' @section Help:
#' For a quick help see the README file:
#' \href{https://raw.githack.com/dppalomar/spectralGraphTopology/master/README.html}{GitHub-README}.
#'
#' For more details see the vignette:
#' \href{https://raw.githack.com/dppalomar/spectralGraphTopology/master/vignettes/SpectralGraphTopology.html}{GitHub-html-vignette}
#'
#' @author Ze Vinicius and Daniel P. Palomar
#'
#' @references
#' S. Kumar, J. Ying, J. V. de Miranda Cardoso, and D. P. Palomar (2019).
#  A unified framework for structured graph learning via spectral constraints.
#' <https://arxiv.org/abs/1904.09792>
#'
#' N., Feiping, W., Xiaoqian, J., Michael I., and H., Heng. (2016).
#' The Constrained Laplacian Rank Algorithm for Graph-based Clustering,
#' AAAI'16. http://dl.acm.org/citation.cfm?id=3016100.3016174
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib spectralGraphTopology
#' @docType package
#' @name spectralGraphTopology-package
NULL
