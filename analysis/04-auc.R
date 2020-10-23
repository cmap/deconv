#
# 
# This script prepares and saves all the figures in the paper
#
#
source("_config.R")

# Colors and labels
COL           <- brewer.pal(7, "Pastel1")  
COL.handle    <- c(brewer.pal(9, "Set1"), "#FF7F00")
plate.labels  <- c(DPK = "compound", LIT = "shRNA")
pert.labels   <- c(BRD = "compound", TRCN = "shRNA")
gene.prop.labels <- c("genes in low proportion", "genes in high proportion")

# Load data sets
speed_acc   <- read.table(path_speed_acc, na.strings="n/a", head = T, sep = "\t")
soln_desc   <- read.table(path_soln_desc, na.strings="n/a", head = T, sep = "\t")
ds_cor      <- parse_gctx(path_ds_cor)
annot       <- read.table(path_annot, head = T, na.strings = "-666", sep = "\t")
gene_annot  <- read.table(path_gene_annot, colClasses = c("gene_id"="character"), head = T)
ds_deconv   <- parse_gctx(path_ds_deconv)

# Process data 
speed_acc <- merge(speed_acc, soln_desc, by="handle")
ds_deconv <- annotate.gct(ds_deconv, annot, dim="col", keyfield="id")
ds_deconv <- annotate.gct(ds_deconv, gene_annot, dim="row", keyfield="gene_id")
names(COL.handle) <- soln_desc$handle

# Correlation

plot.auc <- function(x.mean, x.meth, x.baseline = 0.005) {
    y.scale <- c(0, max(x.mean) + 0.002 - x.baseline)
    b <- barplot(x.mean - x.baseline, axisnames = F
        , col = COL[factor(x.meth)], axes = F, border = "white", ylim = y.scale)
    txt <- paste(1:10, x.meth)
    txt[10] <- "dpeak*"
    text(x = b, y = 0.001, txt, adj = 0, srt = 90)
    values.txt <- sprintf("%.1f", 100*x.mean)
    text(x = b, y = x.mean - x.baseline, values.txt, pos = 3, srt = 0, cex = .75, xpd = T)
}
plot.kd <- function(speed_acc) {
  d <- speed_acc %>% filter(!is.na(kd_precision))
  plot(NA, NA, ylim = c(0.65, 0.80), xlim = c(0.025, 0.04), ann = FALSE, axes = FALSE)
  with(d, abline(v = kd_precision[handle=="benchmark"],   lty = 3))
  with(d, abline(h = kd_success_freq[handle=="benchmark"], lty = 3))
#   rect(x.axis[2], y.axis[2], 1, 1, col = rgb(0,0,1, 0.05), border = "white")
  COL <- c(rgb(0,0,1,0.5), rgb(1,0,0,0.75))
  points(kd_success_freq ~ kd_precision, data = d
        , cex = 2, pch = 21, bg = COL[factor(d$handle=="benchmark")], col = "white")
  prec <- d$kd_precision 
  sf <- d$kd_success_freq
  h <- paste(d$rank, d$method_short)
  h[grep("bench", h)] <- "dpeak*"
  pos <- rep(4, length(h))
  pos[1:3] <- 3
  pos[10] <- 1
  pos[6] <- 1
  text(prec, sf, h, pos = pos, xpd = T)
  axis(1)
  axis(2)
  box(bty = "l")
  title(xlab = "precision (TP/TP+FP)", ylab = "recall (TP/TP+FN)")
}


#+ auc, fig.width = 2.5*3 + 0.5, fig.height = 2.5 
COL           <- brewer.pal(7, "Paired")  
PCH <- rep(c(19, 15, 17, 25, 23),2)
fig.letters <- letters
lettering <- function() {
  title(list(fig.letters[1], cex = 2), adj = 0)
  fig.letters[-1]
}
y.range <- c(0.87, 0.935)#range(speed_acc$auc)
layout(matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4, byrow= T), widths = c(0.30,0.30,0.30,0.1))
  #par(las = 1, mfrow = c(1, 1), mar = c(1,2,1,1))
  for(p in c("DPK", "LIT")) {
    d <- speed_acc %>% 
         filter(grepl(p, plate)) %>% 
         arrange(rank)
    plot(c(1, 10),y.range, type = "n", ann = F, axes = F)
    if (p == "DPK") fig.letters <- lettering()
    abline(h = axTicks(2), lty = 3, col = gray(.74), lwd = 0.5)
    points(1:length(d$auc), d$auc, cex = 1.5
      , pch = PCH[factor(d$method_short)]
      , col = COL[factor(d$method_short)]
      )
    axis(1,  at = c(1:10), labels = c(1:9, "BM"), cex.axis = 0.75)
    axis(2, lwd = 0)
    title(xlab = "final ranking", ylab = "AUC extreme modulations")
    mtext(ifelse(p == "DPK", "compounds", "shRNAs"), side = 3, adj = 1, cex = 0.75)
  }  
  precision <- speed_acc$kd_precision
  recall <- speed_acc$kd_success_freq
  labels <- speed_acc$method_short
  labels_ranks <- speed_acc$rank  
  plot(precision, recall
      , cex = 1.5, axes = F, ann = T
      , pch = PCH[factor(labels)]
      , col = COL[factor(labels)]
      , ylim = c(0.65, 0.8), xlim = c(0.03, 0.038)
      )
  fig.letters <- lettering()
  mtext("gene knockdowns", side = 3, adj = 1, cex = 0.75)
  ok <- which(labels_ranks == "benchmark")
  abline(v = na.omit(precision[ok]), lty = 3, lwd = 0.5)
  abline(h = na.omit(recall[ok]), lty = 3, lwd = 0.5)
  text.pos <- rep(1, length(labels_ranks))
  text.pos[3] <- 2; text.pos[17] <- 3; text.pos[10] <- 3
  text(precision, recall, labels_ranks, pos = text.pos, offset = 1, col = 2)
  box(bty = "l")
  axis(1); axis(2)
  old.par <- par()
  par(mar = c(0,0,0,0))
  plot(c(0,1), c(0,1), type = "n", axes = F, ann = F)
  legend("center", legend = levels(factor(labels)), col = COL, pch = PCH, bty = "n", xpd = T)
  par(old.par)

 