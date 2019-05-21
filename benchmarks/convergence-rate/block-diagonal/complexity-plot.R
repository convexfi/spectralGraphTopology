library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(scales)
library(matrixStats)

time_list <- readRDS(file = "time-complexity.RDS")

n_lists <- length(time_list)
time <- c()
for (i in 1:n_lists) {
  len <- length(time_list[[i]])
  time <- c(time, time_list[[i]][len])
}
n <- c(20, 40, 80, 120, 160, 200)
gr = .5 * (1 + sqrt(5))
lty = c(1, 3)
colors <- c("#ff5252")
setEPS()
cairo_ps("complexity.ps", family = "Times", height = 5,
         width = gr * 3.5, fallback_resolution = 1500, pointsize = 14)
plot(n, time, type = "l", lty = lty[1], col = colors[1],
     xlab = "p", ylab = "CPU time", lwd = 2)
dev.off()
