library(extrafont)
library(latex2exp)

gamma_rel_err <- readRDS("gamma_rel_err.rds")
gamma_fscore <- readRDS("gamma_fscore.rds")
gamma_set <- c(1, 10, 20, 50, 1e2, 1e3, 1e4, 1e5)
gamma_size <- c(1:length(gamma_set))

colors <- c("#843B62", "#F67E7D", "#FFB997")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error_gamma.ps", family = "ComputerModern", height = 5.5, width = gr * 4.5, pointsize = 14)
plot(gamma_size, gamma_rel_err[,1], type = "b", pch=15, cex=.75, ylim=c(0, max(gamma_rel_err)),
     xlab = TeX("$\\gamma$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
grid()
lines(gamma_size, gamma_rel_err[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(gamma_size, gamma_rel_err[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
axis(side = 1, at = gamma_size, labels = gamma_set)
legend("topright", legend = c(TeX("$n / p = 10$"), TeX("$n / p = 100$"), TeX("$n / p = 1000$")),
       col=colors, pch=c(15, 16, 17), lty=c(1, 1, 1), bty="n")
dev.off()
embed_fonts("relative_error_gamma.ps", outfile="relative_error_gamma.ps")
setEPS()
postscript("fscore_gamma.ps", family = "ComputerModern", height = 5.5, width = gr * 4.5, pointsize = 14)
plot(gamma_size, gamma_fscore[,1], ylim=c(min(gamma_fscore), 1.), xlab = TeX("$\\gamma$"),
     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
grid()
lines(gamma_size, gamma_fscore[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(gamma_size, gamma_fscore[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
axis(side = 1, at = gamma_size, labels = gamma_set)
legend("bottomright", legend = c(TeX("$n / p = 10$"), TeX("$n / p = 100$"), TeX("$n / p = 1000$")),
       col=colors, pch=c(15, 16, 17), lty=c(1, 1, 1), bty="n")
dev.off()
embed_fonts("fscore_gamma.ps", outfile="fscore_gamma.ps")
