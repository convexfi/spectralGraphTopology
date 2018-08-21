#' Constructs a block diagonal matrix from a list of square matrices
#' @param ... list of matrices or individual matrices
#' @return block diagonal matrix

blockDiag <- function(...) {
  return(blockDiagCpp(list(...)))
}
