library(pixmap)

n <- 1
S <- matrix(0, 400, 400)
for (i in c(1:400)) {
  print(i)
  for (j in c(i:400)) {
    X <- read.pnm(paste0(i, ".pgm"))@grey
    Y <- read.pnm(paste0(j, ".pgm"))@grey
    S[i, j] <- sum(diag(t(X) %*% Y)) / sqrt(sum(diag(t(X) %*% X)) * sum(diag(t(Y) %*% Y)))
  }
}

S <- S + t(S) - diag(diag(S))

saveRDS(S, "correlation-matrix.rds")
