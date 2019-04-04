library(extrafont)
library(latex2exp)
library(matrixStats)
library(scales)
ratios <- c(5, 10, 30, 100, 250, 500, 1000)
n_ratios <- c(1:length(ratios))

file_err_naive <- readRDS(file = "block-rel-err-naive.rds")
file_err_sgl   <- readRDS(file = "block-rel-err-SGL.rds")
file_err_qp    <- readRDS(file = "block-rel-err-QP.rds")
file_fs_naive  <- readRDS(file = "block-fscore-naive.rds")
file_fs_sgl    <- readRDS(file = "block-fscore-SGL.rds")
file_fs_qp     <- readRDS(file = "block-fscore-QP.rds")

print(file_err_naive)

rel_err_naive <- colMeans(file_err_naive)
rel_err_qp    <- colMeans(file_err_qp)
rel_err_sgl   <- colMeans(file_err_sgl)
fscore_naive  <- colMeans(file_fs_naive)
fscore_qp     <- colMeans(file_fs_qp)
fscore_sgl    <- colMeans(file_fs_sgl)

sd_rel_err_naive <- colSds(file_err_naive)
sd_rel_err_qp    <- colSds(file_err_qp)
sd_rel_err_sgl   <- colSds(file_err_sgl)
sd_fscore_naive  <- colSds(file_fs_naive)
sd_fscore_qp     <- colSds(file_fs_qp)
sd_fscore_sgl    <- colSds(file_fs_sgl)

colors <- c("#0B032D", "#843B62", "#6ABA81")
pch <- c(15, 7, 9)
lty <- c(1, 1, 1)
legend <- c("Naive", "QP", "SGLA (proposed)")
xlab <- TeX("$\\mathit{n} / \\mathit{p}$")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/relative_error_block_bipartite.ps", family = "ComputerModern", height = 5, width = gr * 3.5, pointsize = 14)
plot(n_ratios, rel_err_naive, type = "b", lty = lty[1], pch=pch[1], cex=.75, ylim=c(0, 1),
     xlab = xlab, ylab = "Average Relative Error", col = colors[1], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(rel_err_naive - sd_rel_err_naive, rev(rel_err_naive + sd_rel_err_naive)),
#        col = alpha(colors[1], alpha = .1), border = NA)
grid()
lines(n_ratios, rel_err_qp, type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(rel_err_qp - sd_rel_err_qp, rev(rel_err_qp + sd_rel_err_qp)),
#        col = alpha(colors[2], alpha = .1), border = NA)
lines(n_ratios, rel_err_sgl, type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(rel_err_sgl - sd_rel_err_sgl, rev(rel_err_sgl + sd_rel_err_sgl)),
#        col = alpha(colors[3], alpha = .1), border = NA)
axis(side = 1, at = n_ratios, labels = ratios)
legend("topright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
dev.off()
embed_fonts("../latex/figures/relative_error_block_bipartite.ps", outfile="../latex/figures/relative_error_block_bipartite.ps")

setEPS()
postscript("../latex/figures/fscore_block_bipartite.ps", family = "ComputerModern", height = 5, width = gr * 3.5, pointsize = 14)
plot(n_ratios, fscore_naive, ylim=c(.2, 1.), xlab = xlab, ylab = "Average F-score", type = "b",
     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(fscore_naive - sd_fscore_naive, rev(fscore_naive + sd_fscore_naive)),
#        col = alpha(colors[1], alpha = .1), border = NA)
grid()
lines(n_ratios, fscore_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(fscore_qp - sd_fscore_qp, rev(fscore_qp + sd_fscore_qp)),
#        col = alpha(colors[2], alpha = .1), border = NA)
lines(n_ratios, fscore_sgl, type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(fscore_sgl - sd_fscore_sgl, rev(fscore_sgl + sd_fscore_sgl)),
#        col = alpha(colors[3], alpha = .1), border = NA)
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
dev.off()
embed_fonts("../latex/figures/fscore_block_bipartite.ps", outfile="../latex/figures/fscore_block_bipartite.ps")
