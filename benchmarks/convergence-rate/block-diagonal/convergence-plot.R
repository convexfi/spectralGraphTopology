library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(scales)

relative_error_list <- readRDS(file = "relerr.RDS")
fscore_list <- readRDS(file = "fscore.RDS")
nll_list <- readRDS(file = "nll.RDS")
objfun_list <- readRDS(file = "objfun.RDS")
constr_list <- readRDS(file = "constr.RDS")
beta <- 1e3
# plot convergence trend
gr = .5 * (1 + sqrt(5))
colors <- c("#706FD3", "#FF5252", "#33D9B2")
pch <- c(15, 16, 17)
lty <- c(1, 2, 3)

setEPS()
cairo_ps("constraint_value.ps", family = "Serif", height = 5, width = gr * 3.5)
constr <- constr_list[[1]]
iter <- c(1:length(constr))
plot(iter, 2 * constr / beta, type = "l", lty = 1, col = alpha("black", alpha = 0.1),
     xlab = "Iteration Number", ylab = "Constraint Value")
abline(h = 0, col = "red", lty = 2)
for (i in 2:length(constr_list)) {
  constr <- constr_list[[i]]
  iter <- c(1:length(constr))
  lines(iter, 2 * constr / beta, type = "l", lty = 1, col = alpha("black", alpha = 0.1))
}
dev.off()
embed_fonts("constraint_value.ps", outfile="constraint_value.ps")

setEPS()
cairo_ps("objective_function.ps", family = "Serif", height = 5, width = gr * 3.5)
objfun <- objfun_list[[1]]
iter <- c(1:length(objfun))
plot(iter, objfun, type = "l", lty = 1, col = alpha("black", alpha = 0.1),
     xlab = "Iteration Number", ylab = "Objective Function")
abline(h = 0, col = "red", lty = 2)
for (i in 2:length(objfun_list)) {
  objfun <- objfun_list[[i]]
  iter <- c(1:length(objfun))
  lines(iter, objfun, type = "l", lty = 1, col = alpha("black", alpha = 0.1))
}
dev.off()
embed_fonts("objective_function.ps", outfile="objective_function.ps")

setEPS()
cairo_ps("relative_error.ps", family = "Serif", height = 5, width = gr * 3.5)
relative_error <- relative_error_list[[1]]
iter <- c(1:length(relative_error))
plot(iter, relative_error, type = "l", lty = 1, col = alpha("black", alpha = 0.1),
     xlab = "Iteration Number", ylab = "Relative Error", ylim = c(.1, .8))
for (i in 2:length(relative_error_list)) {
  relative_error <- relative_error_list[[i]]
  iter <- c(1:length(relative_error))
  lines(iter, relative_error, type = "l", lty = 1, col = alpha("black", alpha = 0.1))
}
dev.off()
embed_fonts("relative_error.ps", outfile="relative_error.ps")

setEPS()
cairo_ps("fscore.ps", family = "Serif", height = 5, width = gr * 3.5)
fscore <- fscore_list[[1]]
iter <- c(1:length(fscore))
plot(iter, fscore, type = "l", lty = 1, col = alpha("black", alpha = 0.1),
     xlab = "Iteration Number", ylab = "F-score", ylim = c(.4, 1))
for (i in 2:length(fscore_list)) {
  fscore <- fscore_list[[i]]
  iter <- c(1:length(fscore))
  lines(iter, fscore, type = "l", lty = 1, col = alpha("black", alpha = 0.1))
}
dev.off()
embed_fonts("fscore.ps", outfile="fscore.ps")

