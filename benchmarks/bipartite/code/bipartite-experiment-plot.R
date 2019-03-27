library(extrafont)
library(latex2exp)
library(matrixStats)
library(scales)
ratios <- c(10, 100, 500, 1000, 5000)
n_ratios <- c(1:length(ratios))

file_err_cgl   <- readRDS(file = "rel_err_CGL.rds")
file_err_naive <- readRDS(file = "rel_err_naive.rds")
file_err_sgl   <- readRDS(file = "rel_err_SGL.rds")
file_err_qp    <- readRDS(file = "rel_err_QP.rds")
file_fs_cgl    <- readRDS(file = "fscore_CGL.rds")
file_fs_naive  <- readRDS(file = "fscore_naive.rds")
file_fs_sgl    <- readRDS(file = "fscore_SGL.rds")
file_fs_qp     <- readRDS(file = "fscore_QP.rds")
file_acc_cgl    <- readRDS(file = "accuracy_CGL.rds")
file_acc_naive  <- readRDS(file = "accuracy_naive.rds")
file_acc_sgl    <- readRDS(file = "accuracy_SGL.rds")
file_acc_qp     <- readRDS(file = "accuracy_QP.rds")
file_recall_cgl    <- readRDS(file = "recal_CGL.rds")
file_recall_naive  <- readRDS(file = "recal_naive.rds")
file_recall_sgl    <- readRDS(file = "recal_SGL.rds")
file_recall_qp     <- readRDS(file = "recal_QP.rds")
file_spec_cgl    <- readRDS(file = "specificity_CGL.rds")
file_spec_naive  <- readRDS(file = "specificity_naive.rds")
file_spec_sgl    <- readRDS(file = "specificity_SGL.rds")
file_spec_qp     <- readRDS(file = "specificity_QP.rds")

rel_err_cgl   <- colMeans(file_err_cgl)
rel_err_naive <- colMeans(file_err_naive)
rel_err_qp    <- colMeans(file_err_qp)
rel_err_sgl   <- colMeans(file_err_sgl)
fscore_cgl    <- colMeans(file_fs_cgl)
fscore_naive  <- colMeans(file_fs_naive)
fscore_qp     <- colMeans(file_fs_qp)
fscore_sgl    <- colMeans(file_fs_sgl)
accuracy_cgl    <- colMeans(file_acc_cgl)
accuracy_naive  <- colMeans(file_acc_naive)
accuracy_qp     <- colMeans(file_acc_qp)
accuracy_sgl    <- colMeans(file_acc_sgl)
recall_cgl    <- colMeans(file_recall_cgl)
recall_naive  <- colMeans(file_recall_naive)
recall_qp     <- colMeans(file_recall_qp)
recall_sgl    <- colMeans(file_recall_sgl)
spec_cgl    <- colMeans(file_spec_cgl)
spec_naive  <- colMeans(file_spec_naive)
spec_qp     <- colMeans(file_spec_qp)
spec_sgl    <- colMeans(file_spec_sgl)

sd_rel_err_cgl   <- colSds(file_err_cgl)
sd_rel_err_naive <- colSds(file_err_naive)
sd_rel_err_qp    <- colSds(file_err_qp)
sd_rel_err_sgl   <- colSds(file_err_sgl)
sd_fscore_cgl    <- colSds(file_fs_cgl)
sd_fscore_naive  <- colSds(file_fs_naive)
sd_fscore_qp     <- colSds(file_fs_qp)
sd_fscore_sgl    <- colSds(file_fs_sgl)
sd_accuracy_cgl    <- colSds(file_acc_cgl)
sd_accuracy_naive  <- colSds(file_acc_naive)
sd_accuracy_qp     <- colSds(file_acc_qp)
sd_accuracy_sgl    <- colSds(file_acc_sgl)
sd_recall_cgl    <- colSds(file_recall_cgl)
sd_recall_naive  <- colSds(file_recall_naive)
sd_recall_qp     <- colSds(file_recall_qp)
sd_recall_sgl    <- colSds(file_recall_sgl)
sd_spec_cgl    <- colSds(file_spec_cgl)
sd_spec_naive  <- colSds(file_spec_naive)
sd_spec_qp     <- colSds(file_spec_qp)
sd_spec_sgl    <- colSds(file_spec_sgl)


colors <- c("#0B032D", "#843B62", "#F67E7D", "#6ABA81")
pch <- c(15, 7, 8, 9)
lty <- c(1, 1, 1, 1)
legend <- c("Naive", "QP", "CGL", "SGA (proposed)")
xlab <- TeX("$\\mathit{n} / \\mathit{p}$")
gr = .5 * (1 + sqrt(5))
setEPS()
cairo_ps("../latex/figures/relative_error_bipartite.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(n_ratios, rel_err_naive, type = "b", lty = lty[1], pch=pch[1], cex=.75, ylim=c(0, .6),
     xlab = xlab, ylab = "Average Relative Error", col = colors[1], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(rel_err_naive - sd_rel_err_naive, rev(rel_err_naive + sd_rel_err_naive)),
        col = alpha(colors[1], alpha = .1), border = NA)
grid()
lines(n_ratios, rel_err_qp, type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(rel_err_qp - sd_rel_err_qp, rev(rel_err_qp + sd_rel_err_qp)),
        col = alpha(colors[2], alpha = .1), border = NA)
lines(n_ratios, rel_err_cgl,  type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(rel_err_cgl - sd_rel_err_cgl, rev(rel_err_cgl + sd_rel_err_cgl)),
        col = alpha(colors[3], alpha = .1), border = NA)
lines(n_ratios, rel_err_sgl, type = "b", lty=lty[4], pch=pch[4], cex=.75, col = colors[4], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(rel_err_sgl - sd_rel_err_sgl, rev(rel_err_sgl + sd_rel_err_sgl)),
        col = alpha(colors[4], alpha = .1), border = NA)
axis(side = 1, at = n_ratios, labels = ratios)
legend("topright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
dev.off()
embed_fonts("../latex/figures/relative_error_bipartite.ps", outfile="../latex/figures/relative_error_bipartite.ps")

setEPS()
cairo_ps("../latex/figures/fscore_bipartite.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(n_ratios, fscore_naive, ylim=c(.55, 1.), xlab = xlab, ylab = "Average F-score", type = "b",
     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(fscore_naive - sd_fscore_naive, rev(fscore_naive + sd_fscore_naive)),
        col = alpha(colors[1], alpha = .1), border = NA)
grid()
lines(n_ratios, fscore_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(fscore_qp - sd_fscore_qp, rev(fscore_qp + sd_fscore_qp)),
        col = alpha(colors[2], alpha = .1), border = NA)
lines(n_ratios, fscore_cgl,  type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(fscore_cgl - sd_fscore_cgl, rev(fscore_cgl + sd_fscore_cgl)),
        col = alpha(colors[3], alpha = .1), border = NA)
lines(n_ratios, fscore_sgl, type = "b", lty=lty[4], pch=pch[4], cex=.75, col = colors[4], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(fscore_sgl - sd_fscore_sgl, rev(fscore_sgl + sd_fscore_sgl)),
        col = alpha(colors[4], alpha = .1), border = NA)
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
dev.off()
embed_fonts("../latex/figures/fscore_bipartite.ps", outfile="../latex/figures/fscore_bipartite.ps")


setEPS()
cairo_ps("../latex/figures/recall_bipartite.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(n_ratios, recall_naive, ylim=c(.7, 1.), xlab = xlab, ylab = "Average Recall", type = "b",
     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(recall_naive - sd_recall_naive, rev(recall_naive + sd_recall_naive)),
        col = alpha(colors[1], alpha = .1), border = NA)
grid()
lines(n_ratios, recall_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(recall_qp - sd_recall_qp, rev(recall_qp + sd_recall_qp)),
        col = alpha(colors[2], alpha = .1), border = NA)
lines(n_ratios, recall_cgl,  type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(recall_cgl - sd_recall_cgl, rev(recall_cgl + sd_recall_cgl)),
        col = alpha(colors[3], alpha = .1), border = NA)
lines(n_ratios, recall_sgl, type = "b", lty=lty[4], pch=pch[4], cex=.75, col = colors[4], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(recall_sgl - sd_recall_sgl, rev(recall_sgl + sd_recall_sgl)),
        col = alpha(colors[4], alpha = .1), border = NA)
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
dev.off()
embed_fonts("../latex/figures/recall_bipartite.ps", outfile="../latex/figures/recall_bipartite.ps")


setEPS()
cairo_ps("../latex/figures/accuracy_bipartite.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(n_ratios, accuracy_naive, ylim=c(.6, 1.), xlab = xlab, ylab = "Average Accuracy", type = "b",
     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(accuracy_naive - sd_accuracy_naive, rev(accuracy_naive + sd_accuracy_naive)),
        col = alpha(colors[1], alpha = .1), border = NA)
grid()
lines(n_ratios, accuracy_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(accuracy_qp - sd_accuracy_qp, rev(accuracy_qp + sd_accuracy_qp)),
        col = alpha(colors[2], alpha = .1), border = NA)
lines(n_ratios, accuracy_cgl,  type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(accuracy_cgl - sd_accuracy_cgl, rev(accuracy_cgl + sd_accuracy_cgl)),
        col = alpha(colors[3], alpha = .1), border = NA)
lines(n_ratios, accuracy_sgl, type = "b", lty=lty[4], pch=pch[4], cex=.75, col = colors[4], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(accuracy_sgl - sd_accuracy_sgl, rev(accuracy_sgl + sd_accuracy_sgl)),
        col = alpha(colors[4], alpha = .1), border = NA)
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
dev.off()
embed_fonts("../latex/figures/accuracy_bipartite.ps", outfile="../latex/figures/accuracy_bipartite.ps")


setEPS()
cairo_ps("../latex/figures/spec_bipartite.ps", family = "Serif", height = 5, width = gr * 3.5)
plot(n_ratios, spec_naive, ylim=c(.5, 1.), xlab = xlab, ylab = "Average Specificity", type = "b",
     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(spec_naive - sd_spec_naive, rev(spec_naive + sd_spec_naive)),
        col = alpha(colors[1], alpha = .1), border = NA)
grid()
lines(n_ratios, spec_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(spec_qp - sd_spec_qp, rev(spec_qp + sd_spec_qp)),
        col = alpha(colors[2], alpha = .1), border = NA)
lines(n_ratios, spec_cgl,  type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(spec_cgl - sd_spec_cgl, rev(spec_cgl + sd_spec_cgl)),
        col = alpha(colors[3], alpha = .1), border = NA)
lines(n_ratios, spec_sgl, type = "b", lty=lty[4], pch=pch[4], cex=.75, col = colors[4], xaxt = "n")
polygon(x = c(n_ratios, rev(n_ratios)), y = c(spec_sgl - sd_spec_sgl, rev(spec_sgl + sd_spec_sgl)),
        col = alpha(colors[4], alpha = .1), border = NA)
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
dev.off()
embed_fonts("../latex/figures/spec_bipartite.ps", outfile="../latex/figures/spec_bipartite.ps")
