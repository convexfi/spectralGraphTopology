library(matrixStats)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(scales)

relative_error_list <- readRDS(file = "relerr.RDS")
fscore_list <- readRDS(file = "fscore.RDS")
nll_list <- readRDS(file = "nll.RDS")
objfun_list <- readRDS(file = "objfun.RDS")
constr_list <- readRDS(file = "constr.RDS")
time_list <- readRDS(file = "time.RDS")

## SGL
n_lists <- length(constr_list)
constr <- constr_list[[1]]
min_len <- length(constr)
for (i in 2:n_lists) {
  constr <- constr_list[[i]]
  len <- length(constr)
  if (len < min_len)
    min_len <- len
}
iter <- c(1:min_len)
relerr <- matrix(0, n_lists, min_len)
fscore <- matrix(0, n_lists, min_len)
nll <- matrix(0, n_lists, min_len)
time <- matrix(0, n_lists, min_len)

for (i in 1:n_lists) {
  relerr[i, ] <- c(relative_error_list[[i]])[1:min_len]
  fscore[i, ] <- c(fscore_list[[i]])[1:min_len]
  nll[i, ] <- c(nll_list[[i]])[1:min_len]
  time[i, ] <- c(time_list[[i]])[1:min_len]
}

avg_relerr <- colMeans(relerr)
avg_fscore <- colMeans(fscore)
avg_nll <- colMeans(nll)
avg_time <- colMeans(time)
std_relerr <- colSds(relerr)
std_fscore <- colSds(fscore)
std_nll <- colSds(nll)

legend <- c("SGL (proposed)")
gr = .5 * (1 + sqrt(5))
colors <- rev(c("#ff793f", "#2c2c54"))
setEPS()
cairo_ps("relative_error.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(avg_time, avg_relerr, type = "l", lty = 1, col = colors[1],
     xlab = "Average CPU time", ylab = "Average Relative Error", ylim = c(0, .8), lwd = 2)
polygon(x = c(avg_time, rev(avg_time)), y = c(avg_relerr - std_relerr, rev(avg_relerr + std_relerr)),
        col = alpha(colors[1], alpha = .1), border = NA)
legend("topright", legend=legend, col=colors, lty = c(1, 3), bty="n", lwd = c(2, 2))
dev.off()
embed_fonts("relative_error.ps", outfile="relative_error.ps")

setEPS()
cairo_ps("fscore.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(avg_time, avg_fscore, type = "l", lty = 1, col = colors[1],
     xlab = "Average CPU time", ylab = "Average F - score", ylim = c(0.1, 0.8), lwd = 2)
polygon(x = c(avg_time, rev(avg_time)), y = c(avg_fscore - std_fscore, rev(avg_fscore + std_fscore)),
        col = alpha(colors[1], alpha = .1), border = NA)
legend("bottomright", legend=legend, col=colors, lty = c(1, 3), bty="n", lwd = c(2, 2))
dev.off()
embed_fonts("fscore.ps", outfile="fscore.ps")

setEPS()
cairo_ps("nll.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(avg_time, avg_nll, type = "l", lty = 1, col = colors[1],
     xlab = "Average CPU time", ylab = "Average Negative Loglikelihood", ylim = c(-25, 3), lwd = 2)
legend("topright", legend=legend, col=colors, lty = c(1, 3), bty="n", lwd = c(2, 2))
dev.off()
embed_fonts("nll.ps", outfile="nll.ps")

