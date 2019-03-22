library(matrixStats)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(scales)

relative_error_list_beta1e1 <- readRDS(file = "relerr_beta1e1.RDS")
relative_error_list_beta1e2 <- readRDS(file = "relerr_beta1e2.RDS")
relative_error_list_beta1e3 <- readRDS(file = "relerr_beta1e3.RDS")
relative_error_list_beta1e4 <- readRDS(file = "relerr_beta1e4.RDS")
fscore_list_beta1e1 <- readRDS(file = "fscore_beta1e1.RDS")
fscore_list_beta1e2 <- readRDS(file = "fscore_beta1e2.RDS")
fscore_list_beta1e3 <- readRDS(file = "fscore_beta1e3.RDS")
fscore_list_beta1e4 <- readRDS(file = "fscore_beta1e4.RDS")
nll_list_beta1e1 <- readRDS(file = "nll_beta1e1.RDS")
nll_list_beta1e2 <- readRDS(file = "nll_beta1e2.RDS")
nll_list_beta1e3 <- readRDS(file = "nll_beta1e3.RDS")
nll_list_beta1e4 <- readRDS(file = "nll_beta1e4.RDS")
objfun_list_beta1e1 <- readRDS(file = "objfun_beta1e1.RDS")
objfun_list_beta1e2 <- readRDS(file = "objfun_beta1e2.RDS")
objfun_list_beta1e3 <- readRDS(file = "objfun_beta1e3.RDS")
objfun_list_beta1e4 <- readRDS(file = "objfun_beta1e4.RDS")

n_lists <- 20
min_len_beta1e1 <- min(length(objfun_list_beta1e1[[1]]))
min_len_beta1e2 <- min(length(objfun_list_beta1e2[[1]]))
min_len_beta1e3 <- min(length(objfun_list_beta1e3[[1]]))
min_len_beta1e4 <- min(length(objfun_list_beta1e4[[1]]))
for (i in 2:n_lists) {
  tmp1 <- objfun_list_beta1e1[[i]]
  tmp2 <- objfun_list_beta1e2[[i]]
  tmp3 <- objfun_list_beta1e3[[i]]
  tmp4 <- objfun_list_beta1e4[[i]]
  len_beta1e1 <- length(tmp1)
  len_beta1e2 <- length(tmp2)
  len_beta1e3 <- length(tmp3)
  len_beta1e4 <- length(tmp4)
  if (len_beta1e1 < min_len_beta1e1)
    min_len_beta1e1 <- len_beta1e1
  if (len_beta1e2 < min_len_beta1e2)
    min_len_beta1e2 <- len_beta1e2
  if (len_beta1e3 < min_len_beta1e3)
    min_len_beta1e3 <- len_beta1e3
  if (len_beta1e4 < min_len_beta1e4)
    min_len_beta1e4 <- len_beta1e4
}
print(min_len_beta1e4)

iter_beta1e1 <- c(1:min_len_beta1e1)
iter_beta1e2 <- c(1:min_len_beta1e2)
iter_beta1e3 <- c(1:min_len_beta1e3)
iter_beta1e4 <- c(1:min_len_beta1e4)
relerr_beta1e1 <- matrix(0, n_lists, min_len_beta1e1)
relerr_beta1e2 <- matrix(0, n_lists, min_len_beta1e2)
relerr_beta1e3 <- matrix(0, n_lists, min_len_beta1e3)
relerr_beta1e4 <- matrix(0, n_lists, min_len_beta1e4)
fscore_beta1e1 <- matrix(0, n_lists, min_len_beta1e1)
fscore_beta1e2 <- matrix(0, n_lists, min_len_beta1e2)
fscore_beta1e3 <- matrix(0, n_lists, min_len_beta1e3)
fscore_beta1e4 <- matrix(0, n_lists, min_len_beta1e4)
nll_beta1e1 <- matrix(0, n_lists, min_len_beta1e1)
nll_beta1e2 <- matrix(0, n_lists, min_len_beta1e2)
nll_beta1e3 <- matrix(0, n_lists, min_len_beta1e3)
nll_beta1e4 <- matrix(0, n_lists, min_len_beta1e4)
objfun_beta1e1 <- matrix(0, n_lists, min_len_beta1e1)
objfun_beta1e2 <- matrix(0, n_lists, min_len_beta1e2)
objfun_beta1e3 <- matrix(0, n_lists, min_len_beta1e3)
objfun_beta1e4 <- matrix(0, n_lists, min_len_beta1e4)
#time <- matrix(0, n_lists, min_len)

for (i in 1:n_lists) {
  relerr_beta1e1[i, ] <- c(relative_error_list_beta1e1[[i]])[1:min_len_beta1e1]
  relerr_beta1e2[i, ] <- c(relative_error_list_beta1e2[[i]])[1:min_len_beta1e2]
  relerr_beta1e3[i, ] <- c(relative_error_list_beta1e3[[i]])[1:min_len_beta1e3]
  relerr_beta1e4[i, ] <- c(relative_error_list_beta1e4[[i]])[1:min_len_beta1e4]
  fscore_beta1e1[i, ] <- c(fscore_list_beta1e1[[i]])[1:min_len_beta1e1]
  fscore_beta1e2[i, ] <- c(fscore_list_beta1e2[[i]])[1:min_len_beta1e2]
  fscore_beta1e3[i, ] <- c(fscore_list_beta1e3[[i]])[1:min_len_beta1e3]
  fscore_beta1e4[i, ] <- c(fscore_list_beta1e4[[i]])[1:min_len_beta1e4]
  nll_beta1e1[i, ] <- c(nll_list_beta1e1[[i]])[1:min_len_beta1e1]
  nll_beta1e2[i, ] <- c(nll_list_beta1e2[[i]])[1:min_len_beta1e2]
  nll_beta1e3[i, ] <- c(nll_list_beta1e3[[i]])[1:min_len_beta1e3]
  nll_beta1e4[i, ] <- c(nll_list_beta1e4[[i]])[1:min_len_beta1e4]
  objfun_beta1e1[i, ] <- c(objfun_list_beta1e1[[i]])[1:min_len_beta1e1]
  objfun_beta1e2[i, ] <- c(objfun_list_beta1e2[[i]])[1:min_len_beta1e2]
  objfun_beta1e3[i, ] <- c(objfun_list_beta1e3[[i]])[1:min_len_beta1e3]
  objfun_beta1e4[i, ] <- c(objfun_list_beta1e4[[i]])[1:min_len_beta1e4]
  #time[i, ] <- c(time_list[[i]])[1:min_len]
}

avg_relerr_beta1e1 <- colMeans(relerr_beta1e1)
avg_relerr_beta1e2 <- colMeans(relerr_beta1e2)
avg_relerr_beta1e3 <- colMeans(relerr_beta1e3)
avg_relerr_beta1e4 <- colMeans(relerr_beta1e4)
avg_fscore_beta1e1 <- colMeans(fscore_beta1e1)
avg_fscore_beta1e2 <- colMeans(fscore_beta1e2)
avg_fscore_beta1e3 <- colMeans(fscore_beta1e3)
avg_fscore_beta1e4 <- colMeans(fscore_beta1e4)
avg_nll_beta1e1 <- colMeans(nll_beta1e1)
avg_nll_beta1e2 <- colMeans(nll_beta1e2)
avg_nll_beta1e3 <- colMeans(nll_beta1e3)
avg_nll_beta1e4 <- colMeans(nll_beta1e4)
avg_objfun_beta1e1 <- colMeans(objfun_beta1e1)
avg_objfun_beta1e2 <- colMeans(objfun_beta1e2)
avg_objfun_beta1e3 <- colMeans(objfun_beta1e3)
avg_objfun_beta1e4 <- colMeans(objfun_beta1e4)
#avg_time <- colMeans(time)
#std_relerr_beta1e1 <- colSds(relerr_beta1e1)
#std_relerr_beta1e2 <- colSds(relerr_beta1e2)
#std_relerr_beta1e3 <- colSds(relerr_beta1e3)
#std_fscore_beta1e1 <- colSds(fscore_beta1e1)
#std_fscore_beta1e2 <- colSds(fscore_beta1e2)
#std_fscore_beta1e4 <- colSds(fscore_beta1e4)
#std_nll_beta1e1 <- colSds(nll_beta1e1)
#std_nll_beta1e2 <- colSds(nll_beta1e2)
#std_nll_beta1e3 <- colSds(nll_beta1e3)
#std_objfun_beta1e1 <- colSds(objfun_beta1e1)
#std_objfun_beta1e2 <- colSds(objfun_beta1e2)
#std_objfun_beta1e3 <- colSds(objfun_beta1e3)

legend <- c(TeX("SGL $\\beta = 10$"), TeX("SGL $\\beta = 100$"), TeX("SGL $\\beta = 1000$"), TeX("SGL $\\beta = 10000$"))
gr = .5 * (1 + sqrt(5))
colors <- c("#B53471", "#2f3542", "#ff7f50", "#5352ed")
lty <- c(1, 2, 3, 4)
setEPS()
cairo_ps("relative_error.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(iter_beta1e1, avg_relerr_beta1e1, type = "l", lty = lty[1], col = colors[1],
     xlab = "Iteration Number", ylab = "Average Relative Error", ylim = c(0, 1))
lines(iter_beta1e2, avg_relerr_beta1e2, type = "l", lty = lty[2], col = colors[2])
lines(iter_beta1e3, avg_relerr_beta1e3, type = "l", lty = lty[3], col = colors[3])
lines(iter_beta1e4, avg_relerr_beta1e4, type = "l", lty = lty[4], col = colors[4])
legend("topright", legend=legend, col=colors, lty = lty, bty="n")
dev.off()
embed_fonts("relative_error.ps", outfile="relative_error.ps")

setEPS()
cairo_ps("fscore.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(iter_beta1e1, avg_fscore_beta1e1, type = "l", lty = lty[1], col = colors[1],
     xlab = "Iteration Number", ylab = "Average F-score", ylim = c(0, 1))
lines(iter_beta1e2, avg_fscore_beta1e2, type = "l", lty = lty[2], col = colors[2])
lines(iter_beta1e3, avg_fscore_beta1e3, type = "l", lty = lty[3], col = colors[3])
lines(iter_beta1e4, avg_fscore_beta1e4, type = "l", lty = lty[4], col = colors[4])
legend("bottomright", legend=legend, col=colors, lty = lty, bty="n")
dev.off()
embed_fonts("fscore.ps", outfile="fscore.ps")

setEPS()
cairo_ps("objfun.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(iter_beta1e1, avg_objfun_beta1e1, type = "l", lty = lty[1], col = colors[1],
     xlab = "Iteration Number", ylab = "Average Objective Function")
lines(iter_beta1e2, avg_objfun_beta1e2, type = "l", lty = lty[2], col = colors[2])
lines(iter_beta1e3, avg_objfun_beta1e3, type = "l", lty = lty[3], col = colors[3])
lines(iter_beta1e4, avg_objfun_beta1e4, type = "l", lty = lty[4], col = colors[4])
legend("topright", legend=legend, col=colors, lty = lty, bty="n")
dev.off()
embed_fonts("objfun.ps", outfile="objfun.ps")

setEPS()
cairo_ps("nll.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(iter_beta1e1, avg_nll_beta1e1, type = "l", lty = lty[1], col = colors[1],
     xlab = "Iteration Number", ylab = "Average Negative Loglikelihood", ylim = c(-132, -116))
lines(iter_beta1e2, avg_nll_beta1e2, type = "l", lty = lty[2], col = colors[2])
lines(iter_beta1e3, avg_nll_beta1e3, type = "l", lty = lty[3], col = colors[3])
lines(iter_beta1e4, avg_nll_beta1e4, type = "l", lty = lty[4], col = colors[4])
legend("topright", legend=legend, col=colors, lty = lty, bty="n")
dev.off()
embed_fonts("nll.ps", outfile="nll.ps")
