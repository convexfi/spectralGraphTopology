library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(scales)
library(matrixStats)

times <- readRDS(file = "time-complexity.RDS")
avg_time <- colMeans(times)

n <- c(20, 40, 80, 120, 160, 200, 400)
gr = .5 * (1 + sqrt(5))
lty = c(1, 3)
colors <- c("#ff5252")
setEPS()
cairo_ps("complexity.ps", family = "Times", height = 5,
         width = gr * 3.5, fallback_resolution = 1500, pointsize = 14)
plot(n, avg_time, type = "l", lty = lty[1], col = colors[1],
     xlab = "number of nodes, p", ylab = "Time until convergence [seconds]", lwd = 2)
dev.off()
