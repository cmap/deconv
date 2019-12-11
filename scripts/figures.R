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


# Get deconv data with annotations 
get.ds_deconv <- function() {
  ds <- parse.gctx("data/DECONV_holdout_n8228x976.gctx")
  annot <- fread("data/holdout_sample_annotations.txt")
  gene_annot <- fread("data/barcode_to_gene_map.txt", colClasses = c("gene_id"="character"))
  ds <- annotate.gct(ds, annot, dim="col", keyfield="id")
  ds <- annotate.gct(ds, gene_annot, dim="row", keyfield="gene_id")
  return(ds)
}

## Inter-replicate variance 
compute.inter.rep <- function() {
  ds <- get.ds_deconv()
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
  return(out)
}

plot.inter.ecdf <- function(out) {
  extract.var <- function(out, handle, rows) {
    out[[handle]]$var[rows, ]
  } 
  par(ann = F)
  for (k in c("high", "low")) {
      library(RColorBrewer)
      b.pal <- brewer.pal(5, "Blues")
      o.pal <- brewer.pal(5, "Oranges")
      COL <- c(high = b.pal[4], low = b.pal[2])
      COL.2 <- c(high = o.pal[4], low = o.pal[2])
      rows <- !high # !high 
      if (k == "high") rows <- high 
      x.bench <- out[["benchmark"]]$var[rows, ]
      x <- extract.var(out,  "gardn999", rows)
      x.mean <- exp(apply(x, 1, mean))  # = 1 average variance by genes
      y.mean <- exp(apply(x.bench, 1, mean))
      x.scale <- c(1, 1.1)
      if (k == "high") plot(NA, NA, ylim = c(0, 1), xlim = x.scale, axes = F)
      plot.stepfun(ecdf(x.mean), do.points = F, lwd = 2, col = COL[k], add = T)
      plot.stepfun(ecdf(y.mean), do.points = F, lwd = 2, col = COL.2[k], add = T)
      box(bty = "l")
      axis(2, at = c(0, 1, 0.5))      
      axis(1, at = c(0, 1, 1.1, 2))
      title(xlab = "inter-replicate variability"
          , ylab = "cumulative fraction of perturbagens")      
      txt.1 <- ifelse(k==1, "high bead prop.", "low bead prop")
      text(x = 0.7, y = 0.25, sprintf("%i genes\n%s", length(x), txt.1), adj = 0)
      legend.txt <- c("winner (high)","winner (low)", "benchmark (high)", "benchmark (low)")
      legend("bottomright", lty = 1, lwd = 2, col = c(COL, COL.2), legend = legend.txt, , bty = "n")
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

# [DEPRECTATED] Function to plot AUC figures
plot.auc.deprec <- function(speed_acc) {
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



# AUC figure 
plot.auc <- function(speed_acc, x.baseline=0.005) {
#  par(las = 1, mfrow = c(1, 2), mar = c(1,4,1,1), oma = c(0,0,4,0))
  COL <-   brewer.pal(7, "Pastel1")
  for(p in c("DPK", "LIT")) {
    d <- speed_acc %>% 
         filter(grepl(p, plate)) %>% 
         arrange(rank)
      x.mean <- d$auc
      x.meth <- d$method_short
#       x.baseline <- min(x.mean) - offset.baseline
      y.scale <- c(0, max(x.mean) + 0.002 - x.baseline)
      barplot(x.mean - x.baseline
            , axisnames = F, col = COL[factor(x.meth)], axes = F, border = "white", ylim = y.scale)
      title(ylab = "AUC", line = 1, outer = T)
      txt <- paste(1:10, x.meth)
      txt[10] <- "dpeak*"
      text(x = b, y = 0.001, txt, adj = 0, srt = 90)
      values.txt <- sprintf("%.1f", 100*x.mean)
      text(x = b, y = x.mean - x.baseline, values.txt, pos = 3, srt = 0, cex = .75, xpd = T)
      mtext(side = 3, ifelse(p == "DPK", "shRNAs", "compounds"))
  }
}

plot.kd <- function(speed_acc) {
  d <- speed_acc %>% 
       filter(!is.na(kd_precision))
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
 
# DEPRECATED Function to plot corrlaion 
plot.cor.barplot <- function(ds, handles, COL, y.scale = c(-0.01, 0.05), y.offset = 0.0075, ...) { 
  par(las = 1, mfrow = c(1, 2), mar = c(1, 4, 2, 1), oma = c(0,0,4,0))
  for (plate in c("DPK", "LIT")) {
    x <- extract.cols(ds, handle = handles, plate)
    x.mean <- apply(x@mat, 2, mean)
    x.baseline <- min(x.mean) - 0.015
    barplot(x.mean - x.baseline, axisnames = F, col = COL, axes = F, border = "white", ...)
    abline(h = pretty(x.mean) - x.baseline, lty = 2, col = gray(.75))
    b <- barplot(x.mean - x.baseline, axisnames = F, col = COL, axes = F, add = T, border = "white", ...)
    title(main = ifelse(plate == "DPK", "shRNAs", "compounds"), line = 1)
    title(ylab = "avg. Spearman's rank correlation")
    axis(2, at = pretty(x.mean) - x.baseline, pretty(x.mean))
    box()
    txt <- paste(1:10, soln_desc$method_short)
    txt[10] <- "dpeak*"
    text(x = b, y = 0.001, txt, adj = 0, srt = 90)
  }
}

# ECDF winner vs benchmark 
plot.cor.winner <-  function(d) {
  d$value <- ifelse(d$value < 0, 0, d$value)
  COL <- c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))
  y <- factor(d$plate.prop[d$handle == "benchmark"])
  x <- d$value[d$handle == "benchmark"] 
  x.winner <- d$value[d$handle == "gardn999"] 
  y.winner <- factor(d$plate.prop[d$handle == "gardn999"])
  par(mfrow = c(2, 2), ann = F)
  for (k in 1:nlevels(y)) {
    x.ecdf <- ecdf(x[y == levels(y)[k]])
    x.ecdf.winner <- ecdf(x.winner[y.winner == levels(y.winner)[k]])
    plot.stepfun(x.ecdf, do.points = F, col = COL[1])
    plot.stepfun(x.ecdf.winner, do.points = F, col = COL[2], add = T)
    abline(h = seq(0.2, 0.8, by = 0.2), lty = 3)
    pr <- ifelse(grepl("^0 ", levels(y)[k]), "low bead proportion", "high bead proportion")
    pl <- ifelse(grepl("DPK", levels(y)[k]), "(shRNAs)", "(compounds)")
    title(paste(pr, pl))
    stat.bench <- mean(x[y == levels(y)[k]])
    stat.winner <- mean(x.winner[y.winner == levels(y.winner)[k]])
    legend.1 <- bquote("benchmark"~bar(rho)==.(stat.bench))
    legend.2 <- bquote("winner"~bar(rho)==.(stat.winner))
    legend("topleft", lty =1 , col = COL, legend = as.expression(c(legend.1, legend.2)), bty = "n")
    title(xlab = "genewise Spearman's rank correlation")
    title(ylab = "cumulative fraction of genes")
  }
}

plot.cor.bars <- function(ds_cor, ds_deconv, p, x.baseline = 0
                          , handles = soln_desc$handle, x.meth = soln_desc$method_short) {
    COL <-   brewer.pal(7, "Pastel1")  
    for(k in c(0, 1)) {    
      # compute 
      ds <- subset.gct(ds_cor, rid = ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == k])
      d <- extract.cols(ds, handle = handles, plate = p)
      x.mean <- apply(d@mat, 2, mean)
      x.se <- apply(d@mat, 2, function(x)sqrt(var(x)/length(x)))
      if (is.null(x.baseline)) x.baseline <- min(x.mean)# - 0.025
      x.lo <- x.mean - x.se - x.baseline
      x.hi <- x.mean + x.se - x.baseline
      y.scale <- c(0, max(x.hi) + 0.005)
      b <- barplot(x.mean - x.baseline, axisnames = F, axes = F, ylim = y.scale
            , col = COL[factor(x.meth)], border = "white")
      segments(x0=b, y0 = x.lo, y1 = x.hi)
      mtext(side = 3, ifelse(k == 1, "genes in high bead proportion", "genes in low bead proportion"))
#         mtext(side = 3, line = 1, ifelse(p == "DPK", "shRNAs", "compounds"), font = 4)
      title(ylab = "average genewise correlation", line = 1, outer = T)
      txt <- paste(1:10, x.meth)
      txt[10] <- "dpeak*"
      text(x = b, y = 0.001, txt, adj = 0, srt = 90)
      values.txt <- sprintf("%.1f", 100*x.mean)
      text(x = b, y = x.hi, values.txt, pos = 3, srt = 0, cex = .75, xpd = T)
    }      
}

plot.ecdf.cor.winner <- function(ds_cor, ds_deconv) {
  par(ann = F)
  for(k in c(0, 1)) {
    for(p in c("DPK", "LIT")) {
      ds_cor.sub <- subset.gct(ds_cor, rid = ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == k])
      d <- extract.cols(ds_cor.sub, c("gardn999", "benchmark"), plate = p)
      x <- ifelse(d@mat[,1]<0, 0, d@mat[,1])
      y <- ifelse(d@mat[,2]<0, 0, d@mat[,2])
      library(RColorBrewer)
      b.pal <- brewer.pal(5, "Blues")
      o.pal <- brewer.pal(5, "Oranges")
      COL <- c(DPK = b.pal[4], LIT = b.pal[2])
      COL.2 <- c(DPK = o.pal[4], LIT = o.pal[2])
      if (p=="DPK") plot(NA,NA, ylim = c(0, 1), xlim = c(0, 1), axes = F)
      plot.stepfun(ecdf(x), do.points = F, lwd = 2, col = COL[p], add = T)
      plot.stepfun(ecdf(y), do.points = F, lwd = 2, col = COL.2[p], add = T)
      box(bty = "l")
      axis(2, at = c(0, 1, 0.5))      
      axis(1, at = c(0, 1, 0.5))
      title(xlab = "genewise Spearman's rank correlation", ylab = "cumulative fraction of genes")      
      txt.1 <- ifelse(k==1, "high bead prop.", "low bead prop")
      text(x = 0.7, y = 0.25, sprintf("%i genes\n%s", length(x), txt.1), adj = 0)
      legend.txt <- c("winner (shRNAs)","winner (compounds)", "benchmark (shRNAs)", "benchmark (compounds)")
      legend("topleft", lty = 1, lwd = 2, col = c(COL, COL.2), legend = legend.txt, , bty = "n")
    }
  }
}

# ================================================================== #
# ================================================================== #
# ================================================================== #

main <- function() {
  
  # Get data 
  speed_acc <- fread("data/holdout_score_time.txt", na.strings="n/a")
  soln_desc <- fread("data/solution_descriptions.txt", na.strings="n/a")
  speed_acc <- merge(speed_acc, soln_desc, by="handle")
  ds_cor <- parse.gctx("data/UNI_DUO_gene_spearman_correlations_holdout_n20x976.gctx")
  ds_deconv <- get.ds_deconv()

  # Correlation by low/high bead proportion 
  ds_cor.low <- subset.gct(ds_cor, rid = ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == 0])
  ds_cor.high <- subset.gct(ds_cor, rid = ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == 1])
  extract.cols <- function(x, handle, plate) {
    j <- which(x@cdesc$handle %in% handle &  grepl(plate, x@cdesc$plate))
    x.sub <- subset.gct(x, cid = j)
    subset.gct(x.sub, cid = match(handle, x.sub@cdesc$handle))
  }


  pdf("outcomes/figures/corr-ecdf-winner.pdf", width = 6, height = 6, onefile = T)
  plot.ecdf(ds_cor.winner, ds_deconv)
  dev.off()

  pdf("outcomes/figures/corr-bars.pdf", width = 8, height = 4, onefile = T)
  for (plate in c("DPK", "LIT")) {
    par(las = 1, mfrow = c(1, 2), mar = c(0,0,1,1), oma = c(1,3,0,0))
    plot.cor.bars(ds_cor, ds_deconv, p = plate, x.baseline = ifelse(plate == "DPK", 0.49, 0.63))
  }
  dev.off()

  pdf("outcomes/figures/auc-finale.pdf", width = 8, height = 4)
  par(las = 1, mfrow = c(1, 2), mar = c(0,0,1,1), oma = c(1,3,0,0))
  plot.auc(speed_acc, x.baseline = 0.86)
  dev.off()

  # Data melt 
  gene_annot <- fread("data/barcode_to_gene_map.txt", colClasses = c("gene_id"="character"))
  d <- melt(ds_cor@mat) %>% 
        rename(gene_id = Var1) %>%
        tidyr::separate(Var2, into = c("plate", "handle"), sep = ":") %>% 
        mutate(gene_id = as.character(gene_id)) %>% 
        mutate(value.bc = ifelse(value < 0.05, 0.05, value)) %>% 
        left_join(gene_annot) %>%
        mutate(plate.prop = paste(high_prop, plate)) %>% 
        as_tibble
  pdf("outcomes/figures/cor_winner.pdf", width = 8, height = 8)
  plot.cor.winner(d)
  dev.off()

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
  par(las = 1, mfrow = c(1, 2), mar = c(0,0,1,1), oma = c(1,3,0,0))
  plot.auc(speed_acc)
  dev.off()
  
  pdf("outcomes/figures/kd.pdf", width = 5, height = 5)
  par(mar = c(4,4,1,1))
  plot.kd(speed_acc)
  dev.off()

  out <- compute.inter.rep()
  
  pdf("outcomes/figures/inter-rep.pdf", width = 5, height = 5, onefile = T)
  plot.inter.ecdf(out)
  dev.off()

}

# main()
