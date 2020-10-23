#
# 
# This script prepares and saves all the figures in the paper
#
#

#+ setup, include = FALSE 
library(cmapR)
library(dplyr, warn = F)
library(RColorBrewer)
library(ggsci)
source("_config.R")

# Load data sets -------------------------------------------------------
speed_acc   <- read.table(path_speed_acc, na.strings="n/a", head = T, sep = "\t")
soln_desc   <- read.table(path_soln_desc, na.strings="n/a", head = T, sep = "\t")
ds_cor      <- parse_gctx(path_ds_cor)
annot       <- read.table(path_annot, head = T, na.strings = "-666", sep = "\t")
gene_annot  <- read.table(path_gene_annot, colClasses = c("gene_id"="character"), head = T)
ds_deconv   <- parse_gctx(path_ds_deconv)

# Process data -------------------------------------------------------
speed_acc <- merge(speed_acc, soln_desc, by="handle")
ds_deconv <- annotate.gct(ds_deconv, annot, dim="col", keyfield="id")
ds_deconv <- annotate.gct(ds_deconv, gene_annot, dim="row", keyfield="gene_id")
names(COL.handle) <- soln_desc$handle

# Correlation

## extract median correlation by high vs. low 
get_median_cor <- function() {
  k.levels <- c("low prop." = 0, "high prop." = 1)
  out <- array(dim = c(488, 2), dimnames = list(NULL, names(k.levels)))
  for (k in seq_along(k.levels)) {
    rids <- ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == k.levels[k]]
    ds <- ds_cor@mat[rids, ]
    out[, k] <- apply(ds, 1, median)
  }
  return(out)
}
med_cor <- get_median_cor()


#+ correlation, fig.width = 3.5*3, fig.height = 3*2, dev = "png", dpi = 300
layout(matrix(c(1,2,3,3,4,5,6,7), nrow = 2, ncol = 4, byrow = T))
par(mar=c(4,4,3,1))
plate.labels  <- c(DPK = "compound", LIT = "shRNA")
fig.letters <- letters

# panel 1 --- compound vs. shrna (benchmark)
plot(c(0, 1), c(0, 5), type = "n", axes = F, ann = F)
legend("topleft", lty = 1:2, plate.labels, bty = "n")
box(bty = "l"); axis(1); axis(2)
title(xlab = "benchmark Spearman's correlation", ylab = "density")
for (plate in c("DPK", "LIT")) {
  ds <- ds_cor@mat[, grep(plate, colnames(ds_cor@mat))]
  colnames(ds) <- gsub(".*:", "", colnames(ds))
  x <- ds[, "benchmark"]#apply(ds, 1, median)
  lines(density(x, bw = 0.05, from = 0, to = 1), main = "", lty = ifelse(plate == "DPK", 1, 2))
}
title(list(fig.letters[1], cex = 2), adj = 0); fig.letters <- fig.letters[-1]

# panel 2 --- high vs. low prop (benchmark)
plot(c(0, 1), c(0, 5), type = "n", axes = F, ann = F)
legend("topleft", lty = 1:2, c("low", "high"), bty = "n")
box(bty = "l"); axis(1); axis(2)
title(xlab = "benchmark Spearman's correlation", ylab = "density")
for (p in c(0,1)) {
  rids <- ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == p]
  ds <- ds_cor@mat[rids, ]
  x <- as.vector(ds[, grepl("benchmark", colnames(ds))])
  lines(density(x, bw = 0.05, kernel = "epanechnikov", from = 0, to = 1), main = "", lty = ifelse(p == 0, 1, 2))
}
ks <- ks.test(med_cor[,1], med_cor[,2]) 
legend("left", title = "Kolmogorov-Smirnov", bty = 'n'
  , legend = ifelse(ks$p.value < 0.01, "p < 0.01", paste("p =", round(ks$p.value,3))))
title(list(fig.letters[1], cex = 2), adj = 0); fig.letters <- fig.letters[-1]

# aggregate results 
ds <- ds_cor@mat
colnames(ds) <- gsub(".*:", "", colnames(ds))
ds <- ds[, soln_desc$handle]
est <- apply(ds, 2, mean)
se <- apply(ds, 2, function(x) sd(x)/sqrt(length(x)))    
col.est <- array(COL[factor(soln_desc$method_short)], dimnames = list(names(est)))
rank.est <- array(c(1:(length(est)-1), "BM"), dimnames = list(names(est)))
method.est <- array(soln_desc$method_short, dimnames = list(names(est)))

est <- est - est["benchmark"]
se <- sqrt(se^2 + se["benchmark"]^2)
est <- est[1:9]
se <- se[1:9]

coefplot <- function(est, se, xlab = "", ylab = "", cex.method = 1, y.range = NULL) {
  if(is.null(y.range)) y.range <- range(c(est + 2.1*se, est-2.1*se))
  x.range <- c(1, length(est) + 0.5)
  par(mar = c(5.1, 4, 3, 1))
  plot(x.range, y.range, axes = F, ann = F, type = "n")
  abline(h = axTicks(2), lty = 3, col = gray(.74), lwd = 0.5)
  points(1:length(est), est, pch = 19, col = COL[factor(soln_desc$method_short)])
  segments(x0 = 1:length(est), y0 = est + se, y1 = est - se
        , col = col.est[names(est)], lwd = 3)
  segments(x0=1:length(est), y0 = est + 1.96*se, y1 = est - 1.96*se
        , col = col.est[names(est)], lwd = 1.5)
  lablist <- paste(rank.est[names(est)], method.est[names(est)], sep ="/")
  axis(1,  at = seq_along(est), labels = lablist, cex.axis = 0.75, tick = F, las = 2, line = 0)
#  text(x = seq_along(est), y = par("usr")[3] - 0.0025, labels = lablist, srt = 30, pos = 2, xpd = TRUE)
  axis(2, at = axTicks(2), labels = round(axTicks(2) * 100), lwd = 0, las = 2)
#  lines(x = x.range, y = c(0, 0), lwd = 2)
  box(bty = "l")
  abline(h = 0, lwd = 1.5, lty = 2)
#   text(x = x.range[1] + 1, y = 0, labels = "benchmark", pos = 1, adj = 0)
  #text(1:length(est), min(est), soln_desc$method_short, cex = 1, pos = 1, xpd = T)
#   text(x = (1:length(est)) - 0.3, y = est, labels = method.est[names(est)]
#       , cex = cex.method, xpd = T, srt = 90, col = col.est[names(est)])
  title(xlab = "final ranking / method", line = 4)
  title(ylab = "mean diff. w/ benchmark")
}

coefplot(est, se)
mtext("full sample", side = 3, adj = 1, cex = 0.75)
title(list(fig.letters[1], cex = 2), adj = 0); fig.letters <- fig.letters[-1]

# Stratified 
for (plate in c("DPK", "LIT")) {
  for (k in c(0, 1)) {
    rids <- ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == k]
    ds <- ds_cor@mat[rids, grep(plate, colnames(ds_cor@mat))]
    colnames(ds) <- gsub(".*:", "", colnames(ds))
    ds <- ds[, soln_desc$handle]
    est <- apply(ds, 2, mean)
    se <- apply(ds, 2, function(x) sd(x)/sqrt(length(x))) 
    # difference
    est <- est - est["benchmark"]
    se <- sqrt(se^2 + se["benchmark"]^2)
    est <- est[1:9]
    se <- se[1:9]
    #
    coefplot(est, se, xlab = "final ranking", ylab = "mean improvement"
            , cex.method = 0.0001, y.range = c(-1/2, 1) * 0.08)  
    if (k==0 & plate =="DPK") title(list(fig.letters[1], cex = 2), adj = 0); fig.letters <- fig.letters[-1]
    mtext(paste(plate.labels[plate], "/", ifelse(k==0, "low prop.", "high prop."))
          , side = 3, adj = 1, cex = 0.75)
    }
}


