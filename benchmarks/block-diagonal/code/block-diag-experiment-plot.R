library(extrafont)
library(latex2exp)
library(matrixStats)
library(scales)
ratios <- c(10, 50, 100, 5e2, 1e3)
n_ratios <- c(1:length(ratios))

file_err_naive <- readRDS(file = "rel-err-naive.rds")
file_err_sgl   <- readRDS(file = "rel-err-SGL.rds")
file_err_qp    <- readRDS(file = "rel-err-QP.rds")
file_fs_naive  <- readRDS(file = "fscore-naive.rds")
file_fs_sgl    <- readRDS(file = "fscore-SGL.rds")
file_fs_qp     <- readRDS(file = "fscore-QP.rds")
file_acc_naive  <- readRDS(file = "accuracy-naive.rds")
file_acc_sgl    <- readRDS(file = "accuracy-SGL.rds")
file_acc_qp     <- readRDS(file = "accuracy-QP.rds")
file_recall_naive  <- readRDS(file = "recall-naive.rds")
file_recall_sgl    <- readRDS(file = "recall-SGL.rds")
file_recall_qp     <- readRDS(file = "recall-QP.rds")
file_spec_naive  <- readRDS(file = "specificity-naive.rds")
file_spec_sgl    <- readRDS(file = "specificity-SGL.rds")
file_spec_qp     <- readRDS(file = "specificity-QP.rds")

rel_err_naive <- colMeans(file_err_naive)
rel_err_qp    <- colMeans(file_err_qp)
rel_err_sgl   <- colMeans(file_err_sgl)
fscore_naive  <- colMeans(file_fs_naive)
fscore_qp     <- colMeans(file_fs_qp)
fscore_sgl    <- colMeans(file_fs_sgl)
accuracy_naive  <- colMeans(file_acc_naive)
accuracy_qp     <- colMeans(file_acc_qp)
accuracy_sgl    <- colMeans(file_acc_sgl)
recall_naive  <- colMeans(file_recall_naive)
recall_qp     <- colMeans(file_recall_qp)
recall_sgl    <- colMeans(file_recall_sgl)
spec_naive  <- colMeans(file_spec_naive)
spec_qp     <- colMeans(file_spec_qp)
spec_sgl    <- colMeans(file_spec_sgl)

sd_rel_err_naive <- colSds(file_err_naive)
sd_rel_err_qp    <- colSds(file_err_qp)
sd_rel_err_sgl   <- colSds(file_err_sgl)
sd_fscore_naive  <- colSds(file_fs_naive)
sd_fscore_qp     <- colSds(file_fs_qp)
sd_fscore_sgl    <- colSds(file_fs_sgl)
sd_accuracy_naive  <- colSds(file_acc_naive)
sd_accuracy_qp     <- colSds(file_acc_qp)
sd_accuracy_sgl    <- colSds(file_acc_sgl)
sd_recall_naive  <- colSds(file_recall_naive)
sd_recall_qp     <- colSds(file_recall_qp)
sd_recall_sgl    <- colSds(file_recall_sgl)
sd_spec_naive  <- colSds(file_spec_naive)
sd_spec_qp     <- colSds(file_spec_qp)
sd_spec_sgl    <- colSds(file_spec_sgl)


colors <- c("#0B032D", "#843B62", "#6ABA81")
pch <- c(15, 7, 9)
lty <- c(1, 1, 1)
legend <- c("Naive", "QP", "SGL (proposed)")
xlab <- TeX("$\\mathit{n} / \\mathit{p}$")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/relative_error_blockdiag.ps", family = "ComputerModern", height = 5, width = gr * 3.5, pointsize = 14)
plot(n_ratios, rel_err_naive, type = "b", lty = lty[1], pch=pch[1], cex=.75, ylim=c(0, .4),
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
embed_fonts("../latex/figures/relative_error_blockdiag.ps",
            outfile="../latex/figures/relative_error_blockdiag.ps")

setEPS()
postscript("../latex/figures/fscore_blockdiag.ps", family = "ComputerModern", height = 5, width = gr * 3.5, pointsize = 14)
plot(n_ratios, fscore_naive, ylim=c(.15, 1.), xlab = xlab, ylab = "Average F-score", type = "b",
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
embed_fonts("../latex/figures/fscore_blockdiag.ps",
            outfile="../latex/figures/fscore_blockdiag.ps")


#setEPS()
#cairo_ps("../latex/figures/recall_blockdiag.ps", family = "Serif", height = 5, width = gr * 3.5, pointsize = 14)
#plot(n_ratios, recall_naive, ylim=c(.7, 1.), xlab = xlab, ylab = "Average Recall", type = "b",
#     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(recall_naive - sd_recall_naive, rev(recall_naive + sd_recall_naive)),
#        col = alpha(colors[1], alpha = .1), border = NA)
#grid()
#lines(n_ratios, recall_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(recall_qp - sd_recall_qp, rev(recall_qp + sd_recall_qp)),
#        col = alpha(colors[2], alpha = .1), border = NA)
#lines(n_ratios, recall_sgl, type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(recall_sgl - sd_recall_sgl, rev(recall_sgl + sd_recall_sgl)),
#        col = alpha(colors[3], alpha = .1), border = NA)
#axis(side = 1, at = n_ratios, labels = ratios)
#legend("bottomleft", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
#dev.off()
#embed_fonts("../latex/figures/recall_blockdiag.ps",
#            outfile="../latex/figures/recall_blockdiag.ps")


#setEPS()
#cairo_ps("../latex/figures/accuracy_blockdiag.ps", family = "Serif", height = 5, width = gr * 3.5, pointsize = 14)
#plot(n_ratios, accuracy_naive, ylim=c(.6, 1.), xlab = xlab, ylab = "Average Accuracy", type = "b",
#     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(accuracy_naive - sd_accuracy_naive, rev(accuracy_naive + sd_accuracy_naive)),
#        col = alpha(colors[1], alpha = .1), border = NA)
#grid()
#lines(n_ratios, accuracy_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(accuracy_qp - sd_accuracy_qp, rev(accuracy_qp + sd_accuracy_qp)),
#        col = alpha(colors[2], alpha = .1), border = NA)
#lines(n_ratios, accuracy_sgl, type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(accuracy_sgl - sd_accuracy_sgl, rev(accuracy_sgl + sd_accuracy_sgl)),
#        col = alpha(colors[3], alpha = .1), border = NA)
#axis(side = 1, at = n_ratios, labels = ratios)
#legend("bottomright", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
#dev.off()
#embed_fonts("../latex/figures/accuracy_blockdiag.ps", outfile="../latex/figures/accuracy_blockdiag.ps")


#setEPS()
#cairo_ps("../latex/figures/spec_blockdiag.ps", family = "Serif", height = 5, width = gr * 3.5, pointsize = 14)
#plot(n_ratios, spec_naive, ylim=c(.8, 1.), xlab = xlab, ylab = "Average Specificity", type = "b",
#     pch=pch[1], lty=lty[1], cex=.75, col = colors[1], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(spec_naive - sd_spec_naive, rev(spec_naive + sd_spec_naive)),
#        col = alpha(colors[1], alpha = .1), border = NA)
#grid()
#lines(n_ratios, spec_qp,   type = "b", lty=lty[2], pch=pch[2], cex=.75, col = colors[2], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(spec_qp - sd_spec_qp, rev(spec_qp + sd_spec_qp)),
#        col = alpha(colors[2], alpha = .1), border = NA)
#lines(n_ratios, spec_sgl, type = "b", lty=lty[3], pch=pch[3], cex=.75, col = colors[3], xaxt = "n")
#polygon(x = c(n_ratios, rev(n_ratios)), y = c(spec_sgl - sd_spec_sgl, rev(spec_sgl + sd_spec_sgl)),
#        col = alpha(colors[3], alpha = .1), border = NA)
#axis(side = 1, at = n_ratios, labels = ratios)
#legend("bottomleft", legend=legend, col=colors, pch=pch, lty=lty, bty="n")
#dev.off()
#embed_fonts("../latex/figures/spec_blockdiag.ps", outfile="../latex/figures/spec_blockdiag.ps")
