#' Variability plot
#' ================

#+ setup, include = FALSE
source("_config.R")

#+ data
gene_annot      <- read.table(file.path(path_data, "barcode_to_gene_map.txt")
                              ,  colClasses = c("gene_id"="character"), head = T)
path_annot      <- file.path(path_data, "holdout_sample_annotations.txt")
path_ds_deconv  <- file.path(path_data, "DECONV_holdout_n8228x976.gctx")
annot           <- read.table(path_annot, head = T, sep = "\t")
cv              <- readRDS(file.path(path_data, "robust_cv.rds"))

#+ plots

## params -------------------------------------------------------
col.vec <- rep("gray", 11)
col.vec[which(dimnames(cv)[[3]] == "gardn999")] <- "navy"
col.vec[which(dimnames(cv)[[3]] == "benchmark")] <- "brown"
x <- count(annot, pert_id, pert_iname)
part_names <-array(x$pert_iname, dimnames = list(x$pert_id))

# Compute MAD 
x <- apply(cv, 3, function(x) apply(x, 2, mad))
par(mfrow = c(3,3), mar = rep(1, 4), oma = c(4, 4, 1, 1))
for (k in setdiff(colnames(x), c("ground-truth","benchmark"))) {
  plot(density(x[,"benchmark"], bw = 0.25), bty = "l", ylim = c(0, .8), xlim = c(1, 8), ann = F, axes = F)
  box(bty = "l")
  lines(density(x[,k], bw = 0.25), lty = 2)
  title(main = k)
}
title(xlab = "MAD across 10 replicates (compounds)", ylab = "density", outer = T, line = 1)

# Compute mean across all genes 
cv.arr <- apply(cv, 3, function(x) apply(x, 2, mean))
cv.arr <- cv.arr[order(apply(cv.arr, 1, mean)), ]
cv.arr <- cv.arr[, -which(colnames(cv.arr) == "ground-truth")]

# compute for high and low props 
gene_hp <- intersect(gene_annot$gene_id[gene_annot$high_prop == 1], rownames(cv))
gene_lp <- intersect(gene_annot$gene_id[gene_annot$high_prop == 0], rownames(cv))
cv_hp.mean <- apply(cv, 3, function(x) apply(x[gene_hp, ], 2, median))
cv_lp.mean <- apply(cv, 3, function(x) apply(x[gene_lp, ], 2, median))

# plot(density(log(cv_hp.mean[,"benchmark"])))
# lines(density(log(cv_lp.mean[,"benchmark"])),lty=2)

#+ replicability, fig.width = 3.5* 4, fig.height = 3.5 + 0.25, echo = FALSE 
layout(matrix(c(1,1,1,1, 2:5), nrow = 2, byrow = T), heights = c(0.2, 0.8))

# legend -------------------------------------------------------
old.par <- par()
par(mar=rep(0, 4))
plot(0,1, type = "n", axes = F, ann = F)
legend("center", c("winner", "benchmark"), lty = 1, lwd = 2, col = c("navy", "brown"), bty = "n", ncol = 2)
par(mar = c(10,4,3,1))

# panels ------------------------------------------------------- 
groups <- cut(apply(cv.arr, 1, mean), breaks = quantile(apply(cv.arr, 1, mean)), include = T)
titles <- c("very low variability (Q1)", "low variability (Q2)", "medium variability (Q3)", "high variability (Q4)")
for(g in levels(groups)) {
  ok <- groups == g
  y <- cv.arr[ok, ]
  y.range <- range(y)
  x.range <- c(1, sum(ok))
  plot(x.range, y.range
    , log = "y", type = "n", bty = "l"
    , xlab = "", ylab = "robust CV (mean across all genes)"
    , axes = F
    )
  axis(2)
  for (k in 1:ncol(y)) { 
    lines(y[, k], lwd = 1.5, col = col.vec[k])
  }
  mtext(side = 1, at = 1:sum(ok), las = 2, cex  = .5
      , gsub("\\(.*\\)", "", part_names[rownames(y)]))
  mtext(titles[1], font=1, adj = 0); titles <- titles[-1]
}
