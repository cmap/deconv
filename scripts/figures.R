#
# 
# This script prepares and saves all the figures in the paper
#
#

library(cmapR)
library(dplyr, warn = F)
library(pROC)
library(parallel)

# Utils 
extract.column <- function(ds, handle, plate) {
  j <- which(ds@cdesc$handle == handle & grepl(plate, ds@cdesc$plate))
  ifelse(ds@mat[, j]<0,0,ds@mat[, j])
}
bottomcode <- function(x, value = 0) {ifelse(x<value, value, x)}

# This plots ecdf of genewise corr by solution (compare against benchmark) 
plot.ecdf.corr <- function(ds_cor, handles, method, plate, x.scale = c(0, 1), y.scale = c(0, 1)) {
  n <- length(handles)
  COL <- c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))
  par(mfrow = c(1, n), mar = c(0,0,2,1), oma = c(4,4,3,0))
  for (i in 1:n) { 
    x <- extract.column(ds_cor, handles[i], plate)
    y <- extract.column(ds_cor, "benchmark", plate)
    plot.stepfun(ecdf(x), ann = F, col = COL[1], ylim = y.scale, xlim = x.scale, axes = F, do.points = F)
    plot.stepfun(ecdf(y), col = COL[2], do.points = F, lty = 3, add = T)
    abline(h = c(0.25,0.75,0.5), lty = 3, col = gray(.5))
    axis(1, at =  c(0, 1, 0.5))
    if (i %% n == 1) axis(2, at = c(0, 1, 0.5))
    txt <- sprintf("%i, %s (%s)", i, handles[i], method[i])
    title(main = txt)
    legend("topleft", col = COL, c("competitor", "benchmark"), bty = "n", pch = 19, text.font = 2)
  }
  title(main = ifelse(plate == "DPK", "shRNA", "Compounds"), xlab = "genewise Spearman's rank correlation", ylab = "fraction of genes", outer = T)
}

# This plots scatter of points 
plot.corr <- function(ds_cor, handles, method, plate = NULL, x.scale = c(0.2, 1)) {
  n <- length(handles)
  COL <- c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))
  par(mfrow=c(1, n), mar = c(0,0,1,1), oma = c(4,4,3,0))
  for (i in 1:n) {
    x <- extract.column(ds_cor, handles[i], plate)
    y <- extract.column(ds_cor, "benchmark", plate)
    z <- COL[factor(y > x)]
    plot(y, x, cex = 2, pch = 21, col = "white", bg = z, ylim = c(0, 1), xlim = c(0, 1), axes = F)
    abline(0, 1, lty = 2, lwd =2)
    txt <- paste("%", sprintf("improved = %0.2f", mean(y < x)))
    txt.2 <- paste("%", sprintf("not improved = %0.2f", 1-mean(y < x)))
    legend("topleft", pch = 21, col = COL, pt.bg = COL, bty = "n", legend = c(txt, txt.2))
    abline(h = seq(0.2, 0.8, by = 0.2), lty = 3)
    if(i %% n == 1) axis(2)
    axis(1)
    title(main = sprintf("%i, %s (%s)", i, handles[i], method[i]))
  }
  title(ylab = "competitor's genewise correlation"
      , xlab = "benchmark's genewise rank correlation", outer = T)
  title(main = paste("Top-four performing methods", ifelse(plate == "DPK", "(shRNA)", "(compounds)")), outer = T)
}

## Inter-replicate variance 
plot.inter.rep <- function(bead.proportion.high = FALSE) {
  ds <- parse.gctx("data/DECONV_holdout_n8228x976.gctx")
  annot <- fread("data/holdout_sample_annotations.txt")
  gene_annot <- fread("data/barcode_to_gene_map.txt", colClasses = c("gene_id"="character"))
  ds <- annotate.gct(ds, annot, dim="col", keyfield="id")
  ds <- annotate.gct(ds, gene_annot, dim="row", keyfield="gene_id")
  # Params 
  handles <- unique(ds@cdesc$handle)
  perts <- unique(ds@cdesc$pert_id)
  gene_id <- ds@rdesc$gene_id
  n.perts <- length(perts) # 121
  n.genes <- length(gene_id) # 976
  high <- ds@rdesc$high_prop == 1
  n.high <- sum(high) # 488 
  # Compute variance for each handle 
  # for each perturbagen  
  out <- mclapply(handles, function(h) {
    out.var <- out.mean <- array(NA, c(n.genes, n.perts))
    obs <- numeric(length(perts))
    for (j in 1:length(perts)) {
       cids <- which(ds@cdesc$pert_id==perts[j] & ds@cdesc$handle==h)
       mat <- subset.gct(ds, cid=cids, rid=gene_id)@mat
       out.mean[, j] <-  apply(log(mat), 1, mean)
       out.var[, j] <-  apply(log(mat), 1, var)     
       obs[j] <-  length(cids)    
    }
    list(mean = out.mean, var = out.var, n = obs)
  })
  names(out) <- handles
  extract.var <- function(out, handle, rows) {
    out[[handle]]$var[rows, ]
  } 
  for (k in FALSE:TRUE) {
    par(mfrow = c(1, 4), mar = c(0,0,1,1), oma = c(4,4,3,0))
    y.scale <- c(0, 50)
    x.scale <- c(0, 0.2)
    COL <- c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))
    rows <- !high # !high 
    if (k) rows <- high 
    x.bench <- out[["benchmark"]]$var[rows, ]
    types <- c("gardn999","Ardavel","mkagenius", "Ramzes2")
    meth <- c("DTR", "GMM", "k-means", "CNN")
    for (j in 1:length(types)){
      x <- extract.var(out, types[j], rows)
      x <- apply(x, 1, mean)  # = 1 average variance by genes
      plot.stepfun(ecdf(x), xlim = x.scale, ann = F, axes = F, col = COL[1])
      plot.stepfun(ecdf(apply(x.bench, 1, mean)), col = COL[2], add = T)
      title(main = sprintf("%i, %s (%s)", j, types[j], meth[j]))
      axis(1)
      abline(h = seq(0.2, 0.8, by = 0.2), lty = 3)
      if (j == 1) axis(2)
    }
    title(ylab = "fraction of gene-pert. combination", xlab = "average inter-replicate variance", outer = T)
  }
}

# Function plots runtime improvement
plot.runtime <- function(speed_acc) {
  speed_acc <- speed_acc %>% 
               group_by(handle, method_short, language, rank) %>% 
               summarize(sec = mean(sec)) %>%
               arrange(rank)
  ok <- !is.na(speed_acc$handle) #speed_acc$handle != "benchmark"  
  bench <- speed_acc$handle == "benchmark"  
  sec <- speed_acc$sec[ok]
  sec.bench <- speed_acc$sec[bench]
  speedup <- sec.bench/sec
  handle <- speed_acc$handle[ok]
  method <- speed_acc$method_short[ok]
  language <- speed_acc$language[ok]
  handle_f <- factor(handle, levels = handle) #, levels = names(sort(sec.mean)))
  # plot params 
  x.scale <- c(2, 700)
  y.scale <- c(1, nlevels(handle_f))
  COL <- c(rgb(1,0,0,0.75), rgb(0,0,1,0.75), rgb(0,1,0,0.5))
  # 
  par(mar = c(1, 12, 4, 1))
  plot(NA, NA, axes = F, ann = F, xlim = x.scale, ylim = y.scale, log = "x")
  abline(v = sec.bench, col = COL[2], lty = 2)
  points(sec, as.numeric(handle_f), col = "white", bg = COL[factor(handle_f == "benchmark")], cex = 2, pch = 21)
  axis(3, at = c(5, 20, 100, sec.bench))
  txt <- sprintf("%s, %s (%s)", 1:nlevels(handle_f), levels(handle_f), method)
  axis(2, at = 1:nlevels(handle_f), txt, las = 2)
  box()
  text(sec, as.numeric(handle_f), sprintf("%0.0fx (%s)", speedup, language), pos = 4)
  title(main = "Runtime per plate (in seconds)", line = 3)
#   title(ylab = "Top nine approaches", line = 8)
}

# Function to plot AUC figures
plot.auc <- function(speed_acc) {
  ds_gardn999_auc <- parse.gctx("data/DPK.CP003_PC3_24H_X1_B42_DE_gardn999.gct")
  ds_uni_auc <- parse.gctx("data/DPK.CP003_PC3_24H_X1_B42_DE_UNI.gct")
  ds_bench_auc <- parse.gctx("data/DPK.CP003_PC3_24H_X1_B42_DE_benchmark.gct")
  auc_df <- speed_acc %>% arrange(plate, rank)
  x.labels <- auc_df$handle[1:10]
  x <- matrix(auc_df$auc, 10, 2)
  par(mar = c(6,4,3,1))
  COL <- c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))
  b <- barplot(x-0.85, beside = T, axes = F, col = rep(COL, each = 10), border = "white")
  mtext(side = 3, line = 1, at = apply(b, 2, mean), c("shRNA", "compounds"))
  axis(2, at = range(x-0.85), round(range(x), 2))
  axis(1, at = b[,1], x.labels, las = 2)
  axis(1, at = b[,2], x.labels, las = 2)
  title(ylab = "AUC")
}

plot.kd <- function(speed_acc) {
  d <- speed_acc %>% filter(!is.na(kd_precision))
  plot(NA, NA, ylim = c(0.68, 0.8), xlim = c(0.028, 0.038), ann = FALSE, axes = FALSE)
  x.axis <- with(na.omit(speed_acc), c(max(kd_precision), kd_precision[handle=="benchmark"]))
  y.axis <- with(na.omit(speed_acc), c(max(kd_success_freq), kd_success_freq[handle=="benchmark"]))
  axis(1, at = x.axis, format(x.axis, digits = 2))
  axis(2, at = y.axis, format(y.axis, digits = 2))
  box()
  rect(x.axis[2], y.axis[2], 1, 1, col = rgb(0,0,1, 0.05), border = "white")
  points(kd_success_freq ~ kd_precision, data = speed_acc, cex = 2, pch = 21, bg = rgb(1,0,0,0.5), col = "white")
  prec <- d$kd_precision 
  sf <- d$kd_success_freq
  h <- d$handle 
  pos <- rep(4, length(h))
  pos[1:3] <- 3
  text(prec, sf, h, pos = pos)
  with(speed_acc, abline(v = kd_precision[handle=="benchmark"], col = gray(.75), lty = 3))
  with(speed_acc, abline(h = kd_success_freq[handle=="benchmark"], col = gray(.75), lty = 3))
#  title(main = "D. Knockdown detection", xlab = "precision (TP/TP+FP)", ylab = "recall (TP/TP+FN)")
}


# ======================================== #
# ======================================== #
# ======================================== #

main <- function() {
  
  # Get data 
  speed_acc <- fread("data/holdout_score_time.txt", na.strings="n/a")
  soln_desc <- fread("data/solution_descriptions.txt", na.strings="n/a")
  speed_acc <- merge(speed_acc, soln_desc, by="handle")
  ds_cor <- parse.gctx("data/UNI_DUO_gene_spearman_correlations_holdout_n20x976.gctx")

  # ECDF of correlation 
  n <- 4
  handles <- head(soln_desc$handle, n)
  method <- head(soln_desc$method_short, n)
  pdf("outcomes/figures/ecdf.corr.pdf", width = 2*n, height = 3, onefile = T)
  for (plate in c("DPK", "LIT")) 
    plot.ecdf.corr(ds_cor, handles, method, plate)
  dev.off()

  # Scatter plot 
  pdf("outcomes/figures/corr.pdf", width = 2*n, height = 3, onefile = T)
  for (plate in c("DPK", "LIT")) 
    plot.corr(ds_cor, handles, method, plate)
  dev.off()

  pdf("outcomes/figures/speed.pdf", width = 8, height = 5)
  plot.runtime(speed_acc)
  dev.off()
  
  pdf("outcomes/figures/auc.pdf", width = 5, height = 5)
  plot.auc(speed_acc)
  dev.off()
  
  pdf("outcomes/figures/kd.pdf", width = 5, height = 5)
  plot.kd(speed_acc)
  dev.off()

  pdf("outcomes/figures/inter-rep.pdf", width = 2.5*4, height = 3, onefile = T)
  plot.inter.rep()
  dev.off()

}

main()
