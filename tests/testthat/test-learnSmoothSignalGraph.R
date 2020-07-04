context("Smooth signal graph estimation")
library(testthat)
library(patrick)
library(spectralGraphTopology)
library(igraph)

set.seed(42)

with_parameters_test_that("we can recover a simple grid graph with Kalofolias method", {
    p <- 16
    grid <- make_lattice(length = sqrt(p), dim = 2)
    E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(grid))
    diag(Ltrue) <- p * diag(Ltrue) / sum(diag(Ltrue))
    diag_Ltrue <- diag(Ltrue)
    row_sums_L <- rowSums(diag(diag(Ltrue)) - Ltrue)
    Ltrue <- diag_Ltrue * Ltrue / row_sums_L
    diag(Ltrue) <- diag_Ltrue
    Wtrue <- diag(diag(Ltrue)) - Ltrue
    X <- t(MASS::mvrnorm(as.integer(100 * p), mu = rep(0, p), Sigma = MASS::ginv(Ltrue)))
    res <- func(X = X)
    expect_true(fscore(Ltrue, res$laplacian, 1e-4) > .7)
  },
  cases(list(func = learn_smooth_graph))
        #,list(func = learn_graph_sigrep))
)
