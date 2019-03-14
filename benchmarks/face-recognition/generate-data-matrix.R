library(pixmap)
library(spectralGraphTopology)

img <- read.pnm(paste0(1, ".pgm"))@grey
n <- dim(img)[1]
m <- dim(img)[2]
Y <- matrix(0, 400, n * m)
Y[1, ] <- c(vec(img))
for (i in c(2:400)) {
  print(i)
  img <- read.pnm(paste0(i, ".pgm"))@grey
  Y[i, ] <- c(vec(img))
}

saveRDS(Y, "data-matrix.rds")
