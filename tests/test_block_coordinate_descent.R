library(testthat)
library(spectralGraphTopology)

test_that("test_U_update_consistency", {
  w <- c(1, 2, 3, 4, 5, 6)
  U <- U_update(w, 4)
  expect_that(all.equal(t(U) %*% U, diag(array(1., 4)),
                         check.attributes = FALSE), is_true())
})
