library(extrafont)
library(latex2exp)

k_set <- c(1:14)
re <- readRDS(file = "relerror.rds")
fs <- readRDS(file = "fscore.rds")
bic <- readRDS(file = "bic.rds")

gr = .5 * (1 + sqrt(5))
colors <- c("#2c3e50", "#e74c3c")
pch <- c(15, 16)
lty <- c(1, 2)
setEPS()
postscript("../latex/figures/performance_model_mismatch.ps", family = "ComputerModern", height = 5, width = gr * 4)
par(mar = c(5, 5, 3, 5))
plot(k_set, fs, type = "b", pch = pch[1], cex=.75, ylim=c(min(fs), max(fs)), lty = lty[1],
     xlab = TeX("$\\mathit{k}$"), ylab = "Average F-score", col = colors[1])
grid()
par(new = TRUE)
plot(k_set, re, type = "b", pch = pch[2], cex=.75, xlab = "", ylab = "", xaxt = "n", lty = lty[2],
     yaxt = "n", col = colors[2])
axis(side = 4)
mtext("Relative Error", side = 4, line = 3, family = "ComputerModern")
legend("topleft", legend = c("F-score", "Relative Error"), col=colors, pch = pch, lty = lty, bty="n")
dev.off()
embed_fonts("../latex/figures/performance_model_mismatch.ps", outfile="../latex/figures/performance_model_mismatch.ps")

setEPS()
postscript("../latex/figures/bic_model_mismatch.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(k_set, bic, type = "b", pch = pch[1], cex=.75, ylim=c(min(bic), max(bic)), lty = lty[1],
     xlab = TeX("$\\mathit{k}$"), ylab = "Average Bayesian Info Criterion", col = colors[1])
dev.off()
embed_fonts("../latex/figures/bic_model_mismatch.ps", outfile="../latex/figures/bic_model_mismatch.ps")
