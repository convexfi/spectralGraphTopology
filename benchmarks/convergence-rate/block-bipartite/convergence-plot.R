library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(scales)
library(matrixStats)

relative_error_list <- readRDS(file = "relerr.RDS")
fscore_list <- readRDS(file = "fscore.RDS")
nll_list <- readRDS(file = "nll.RDS")
objfun_list <- readRDS(file = "objfun.RDS")
time_list <- readRDS(file = "time.RDS")

n_lists <- length(objfun_list)
objfun <- objfun_list[[1]]
min_len <- length(objfun)
for (i in 2:n_lists) {
  objfun <- objfun_list[[i]]
  len <- length(objfun)
  if (len < min_len)
    min_len <- len
}
iter <- c(1:min_len)
relerr <- matrix(0, n_lists, min_len)
fscore <- matrix(0, n_lists, min_len)
nll <- matrix(0, n_lists, min_len)
objfun <- matrix(0, n_lists, min_len)
time <- matrix(0, n_lists, min_len)

for (i in 1:n_lists) {
  relerr[i, ] <- c(relative_error_list[[i]])[1:min_len]
  fscore[i, ] <- c(fscore_list[[i]])[1:min_len]
  nll[i, ] <- c(nll_list[[i]])[1:min_len]
  objfun[i, ] <- c(objfun_list[[i]])[1:min_len]
  time[i, ] <- c(time_list[[i]])[1:min_len]
}

avg_relerr <- colMeans(relerr)
avg_fscore <- colMeans(fscore)
avg_nll <- colMeans(nll)
avg_objfun <- colMeans(objfun)
avg_time <- colMeans(time)
std_relerr <- colSds(relerr)
std_fscore <- colSds(fscore)
std_nll <- colSds(nll)
std_objfun <- colSds(objfun)

legend <- c("Relative Error", "F-score")
gr = .5 * (1 + sqrt(5))
lty = c(1, 3)
colors <- c("#ff5252", "black")
setEPS()
postscript("convergence-block-bipartite.ps", family = "ComputerModern", height = 5,
         width = gr * 3.5, pointsize = 14)
par(mar = c(5, 5, 3, 5))
plot(iter, avg_relerr, type = "l", lty = lty[1], col = colors[1],
     xlab = "Iteration number", ylab = "Relative Error", ylim = c(.04, .2), lwd = 2)
#polygon(x = c(iter, rev(iter)), y = c(avg_relerr - std_relerr, rev(avg_relerr + std_relerr)),
#        col = alpha(colors[1], alpha = .3), border = NA)
par(new = TRUE)
plot(iter, avg_fscore, type = "l", lty = lty[2], col = colors[2],
     xlab = "", ylab = "", ylim = c(0.45, .95), lwd = 2, xaxt = "n", yaxt = "n")
#polygon(x = c(iter, rev(iter)), y = c(avg_fscore - std_fscore, rev(avg_fscore + std_fscore)),
#        col = alpha(colors[2], alpha = .3), border = NA)
axis(side = 4)
mtext("F-score", side = 4, line = 3)
legend("right", legend=rev(legend), col=rev(colors), lty = rev(lty), bty="n", lwd = c(2, 2))
dev.off()
embed_fonts("convergence-block-bipartite.ps", outfile="convergence-block-bipartite.ps")

legend <- c("Avg. Negative Loglikelihood", "Avg. Objective Function")
setEPS()
cairo_ps("objective-trend-block-bipartite.ps", family = "ComputerModern", height = 5,
         width = gr * 3.5, pointsize = 14)
par(mar = c(5, 5, 3, 5))
plot(iter, avg_nll, type = "l", lty = lty[1], col = colors[1],
     xlab = "Iteration number", ylab = "Avg. Negative Loglikelihood", lwd = 2)
#polygon(x = c(iter, rev(iter)), y = c(avg_nll - .5*std_nll, rev(avg_nll + .5*std_nll)),
#        col = alpha(colors[1], alpha = .3), border = NA)
par(new = TRUE)
plot(iter, avg_objfun, type = "l", lty = lty[2], col = colors[2],
     xlab = "", ylab = "", lwd = 2, xaxt = "n", yaxt = "n")
#polygon(x = c(iter, rev(iter)), y = c(avg_objfun - std_objfun, rev(avg_objfun + std_objfun)),
#        col = alpha(colors[2], alpha = .3), border = NA)
axis(side = 4)
mtext("Avg. Objective Function", side = 4, line = 3)
legend("topright", legend=legend, col=colors, lty = lty, bty="n", lwd = c(2, 2))
dev.off()
