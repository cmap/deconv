#
# 
# This script prepares and saves all the figures in the paper
#
#
usage <- function() {
  message("Rscript figures.R [path to data]")
}
options(warn = -1)
library(cmapR)
library(dplyr, warn = F)
library(parallel)
library(RColorBrewer)

# Colors and labels # ======================================== #
COL           <- brewer.pal(7, "Pastel1")  
COL.handle    <- c(brewer.pal(9, "Set1"), "#FF7F00")
plate.labels  <- c(DPK = "compound", LIT = "shRNA")
pert.labels   <- c(BRD = "compound", TRCN = "shRNA")
gene.prop.labels <- c("genes in low proportion", "genes in high proportion")
wd <- getwd()

# Functions # ======================================== #
bottomcode <- function(x, value = 0) {
  ifelse(x<value, value, x)
}
std.err <- function(x) {
  sqrt(var(as.vector(x)) / length(as.vector(x)))
}
extract.cols <- function(x, handle, plate) {
  j <- which(x@cdesc$handle %in% handle &  grepl(plate, x@cdesc$plate))
  x.sub <- subset.gct(x, cid = j)
  subset.gct(x.sub, cid = match(handle, x.sub@cdesc$handle))
}
pval2stars <- function(x) {
  ifelse(x < 0.001, "***", ifelse(x < 0.05, "**", ifelse(x < 0.1, "*", "")))
}


cor.coefplot <- function(d, x.scale = NULL) {
  ols <- lm(d@mat ~ 1)
  est <- sapply(coef(summary(ols)), function(x) x[1])
  est.se <- sapply(coef(summary(ols)), function(x) x[2])
  est.ci.long <- cbind(est + 2*est.se, est - 2*est.se)
  est.ci <- cbind(est + est.se, est - est.se)
  x.ks <- apply(d@mat, 2, wilcox.test, y = d@mat[,10])
  if (is.null(x.scale)) x.scale <- range(est.ci.long)
  par(ann = F, mar = c(4,4,3,1))
  plot(NA,NA, xlim = x.scale, ylim = c(1,10), axes = F)
  abline(h=1:10, lty = 3, col = gray(.9))
  points(est, 1:length(est), pch = 19, col = COL[factor(x.meth)])
#   segments(x0 = est.ci[,1], x1 = est.ci[,2], y0 = 1:length(est)
#     , col = COL[factor(x.meth)], lwd = 2)
  segments(x0 = est.ci.long[,1], x1 = est.ci.long[,2], y0 = 1:length(est)
    , col = COL[factor(x.meth)])
#   arrows(x0 = est.ci.long[,1], x1 = est.ci.long[,2], y0 = 1:length(est)
#     , code = 3, angle = 90, length = .1, col = COL[factor(x.meth)])
  sig.stars <- pval2stars(sapply(x.ks, function(x) x$p.value)) 
  txt.short <- x.meth #paste(1:10, x.meth) 
#   txt.short[10] <- "dpeak"
  txt <- sprintf("%s %0.3f %s", txt.short, est, sig.stars)
  text(x = x.scale[1], y = 1:length(est), txt, adj = c(0,0)) 
  #, xpd = T, col = COL[factor(x.meth)])
  axis(2, at = c(1:10), las = 1, c(1:9, "dpeak"))
  axis(1)
  box(bty = "l")
}
      cor.coefplot(d, x.scale = c(0, 0.15) + ifelse(plate=="DPK", 0.47, 0.62))


# Plot correlation bars 
plot.cor.bars <- function(ds, p, x.baseline = 0, handles = soln_desc$handle
                    , x.meth = soln_desc$method_short) {
    d <- extract.cols(ds, handle = handles, plate = p)
    d_bench <- extract.cols(ds, handle = "benchmark", plate = p)
    x.mean <- apply(d@mat, 2, mean)
    x.se <- apply(d@mat, 2, function(x) sqrt(var(x)/length(x)))
    x.ks <- apply(d@mat, 2, wilcox.test, y = d_bench@mat)
    if (is.null(x.baseline)) x.baseline <- min(x.mean)# - 0.025
    x.lo <- x.mean - x.se - x.baseline
    x.hi <- x.mean + x.se - x.baseline
    y.scale <- c(0, .75 - x.baseline) #c(0, max(x.hi) + 0.005)
    b <- barplot(x.mean - x.baseline, axisnames = F, axes = F, ylim = y.scale
          , col = COL[factor(x.meth)]
          , border = "white")
    arrows(x0 = b, y0 = x.lo, y1 = x.hi, code = 3, angle = 90, length = 0.1)
    y.ticks <- seq(par("usr")[3], par("usr")[4], length = 2)
    axis(2, at = y.ticks, round(y.ticks + x.baseline, 2))
    txt <- paste(1:10, x.meth)
    txt[10] <- "dpeak"
    sig.stars <- pval2stars(sapply(x.ks, function(x) x$p.value))
    values.txt <- sprintf("%.3f %s", x.mean, sig.stars)
    text(x = b, y = 0.0015, txt, adj = 0, srt = 90)
    text(x = b, y = 0.1, values.txt, adj = 0, srt = 90)
}
#       plot.cor.bars(ds, p = plate, x.baseline = ifelse(plate == "DPK", 0.49, 0.63)
#                   , handles = soln_desc$handle, x.meth = soln_desc$method_short)


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
compute.inter.var <- function(ds) {
  handles <- unique(ds@cdesc$handle)
  perts <- unique(ds@cdesc$pert_id)
  plates <- unique(ds@cdesc$det_plate)
  gene_id <- ds@rdesc$gene_id
  out.var <- mclapply(handles, function(h) {
    out <- array(dim = c(length(gene_id), length(perts))
          , dimnames = list(gene_id, perts))
     for (j in 1:length(perts)) {      
        cids <- which(ds@cdesc$pert_id==perts[j] & ds@cdesc$handle==h)
        mat <- subset.gct(ds, cid = cids, rid = gene_id)@mat
        out[gene_id, j] <-  apply(mat, 1, IQR) / apply(mat, 1, median)   # coefficient of variation        
     }
     return(out)
  })
  names(out.var) <- handles
  return(out.var)
}

plot.inter.ecdf <- function(out) {
  extract.var <- function(out, handle, rows) {
    out[[handle]]$var[rows, ]
  } 
  par(ann = F)
  mean.arr <- array(NA, dim = c(11, 2), dimnames = list(handles, c("high","low")))
  se.arr <- mean.arr
  for (handle in handles) 
  for (k in c("high", "low")) {
      library(RColorBrewer)
      b.pal <- brewer.pal(5, "Blues")
      o.pal <- brewer.pal(5, "Oranges")
      COL <- c(high = b.pal[4], low = b.pal[2])
      COL.2 <- c(high = o.pal[4], low = o.pal[2])
      rows <- !high # !high 
      if (k == "high") rows <- high 
      x.bench <- out[["benchmark"]]$var
      col.indices <- 1:ncol(x.bench)
#      x.bench <- out[["benchmark"]]$var[rows, col.indices] # more than 4
      x <- extract.var(out,  handle, rows)[, col.indices]
      x.mean <- exp(apply(x, 1, mean))  # = 1 average variance by genes
      y.mean <- exp(apply(x.bench, 1, mean))
      mean.arr[handle, k] <- mean(x.mean)
      se.arr[handle, k] <- sqrt(var(x.mean) / length(x.mean))
      x.scale <- c(1, 1.1)
      plot(NA, NA, ylim = c(0, 1), xlim = x.scale, axes = F)
      plot.stepfun(ecdf(x.mean), do.points = F, lwd = 2, col = COL[k], add = T)
      plot.stepfun(ecdf(y.mean), do.points = F, lwd = 2, col = COL.2[k], add = T)
      box(bty = "l")
      axis(2, at = c(0, 1, 0.5))      
      axis(1, at = c(0, 1, 1.1, 2))
      title(xlab = "inter-replicate variability"
          , ylab = "cumulative fraction of perturbagens")      
      txt.1 <- ifelse(k==1, "high bead prop.", "low bead prop")
      text(x = 0.7, y = 0.25, sprintf("%i genes\n%s", length(x), txt.1), adj = 0)
      legend.txt <- c(paste(handle, c("(high)", "(low)")), "benchmark (high)", "benchmark (low)")
      legend("bottomright", lty = 1, lwd = 2, col = c(COL, COL.2), legend = legend.txt, , bty = "n")
    }
  return(list(mean = mean.arr, se = se.arr))
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

plot.winner <- function(data) {
  COL <- list(brewer.pal(3, "Blues"), brewer.pal(3, "Oranges"))
  plates <- c("DPK","LIT")
  handles <- c("gardn999",  "benchmark")
  meth <- c("winner", "dpeak*")
  for (k in c(0,1)) {
    h.labels <- rep()
    plot(NA, NA, xlim = c(0, 1), ylim =c(0,1), axes=F, ann =F)
    for (p in 1:2) {
      for(h in 1:length(handles)) {
        x <- filter(data, handle == handles[h], high_prop == k, grepl(plates[p], plate))$value
        plot.stepfun(ecdf(ifelse(x<0,0,x)), do.points = F, add = T, col = COL[[p]][h+1], lwd = 2)
        axis(2, at = c(0, 1, 0.5))      
        axis(1, at = c(0, 1, 0.5))
        title(xlab = "genewise Spearman's rank correlation", ylab = "cumulative fraction of genes")      
        box(bty="l")
        h.labels <- c(h.labels, paste(meth[h], ifelse(plates[p]=="LIT","(shRNA)","(compound)")))
      }
      mtext(side = 3, ifelse(k == 1, "genes in high bead proportion", "genes in low bead proportion"))
      h.colors <- c(COL[[1]][2:3],COL[[2]][2:3])
      legend("topleft",lty=1,lwd=2,col=h.colors,h.labels,bty="n")
    }
  }
}

plot.ecdf.winner.deprec <- function(ds_cor, ds_deconv) {
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

#- Main function
main <- function() {
  args <- commandArgs(trailingOnly=TRUE) 
  path_out <- args[1]
  path_data <- ifelse(is.na(args[2]), "data", args[2])
  
  path_speed_acc <- file.path(path_data, "holdout_score_time.txt")
  path_soln_desc <- file.path(path_data, "solution_descriptions.txt")
  path_annot <- file.path(path_data, "holdout_sample_annotations.txt")
  path_gene_annot <- file.path(path_data, "barcode_to_gene_map.txt")
  path_ds_cor <- file.path(path_data, "UNI_DUO_gene_spearman_correlations_holdout_n20x976.gctx")
  path_ds_deconv <- file.path(path_data, "DECONV_holdout_n8228x976.gctx")

  # Load data sets
  speed_acc <- fread(path_speed_acc, na.strings="n/a")
  soln_desc <- fread(path_soln_desc, na.strings="n/a")
  ds_cor <- parse.gctx(path_ds_cor)
  annot <- fread(path_annot)
  gene_annot <- fread(path_gene_annot, colClasses = c("gene_id"="character"))
  ds_deconv <- parse.gctx(path_ds_deconv)

  # Process data 
  speed_acc <- merge(speed_acc, soln_desc, by="handle")
  ds_deconv <- annotate.gct(ds_deconv, annot, dim="col", keyfield="id")
  ds_deconv <- annotate.gct(ds_deconv, gene_annot, dim="row", keyfield="gene_id")
  names(COL.handle) <- soln_desc$handle
  cor_df <- melt(ds_cor@mat) %>% 
            rename(gene_id = Var1) %>%
            tidyr::separate(Var2, into = c("plate", "handle"), sep = ":") %>% 
            mutate(gene_id = as.character(gene_id)) %>% 
            mutate(value.bc = ifelse(value < 0.05, 0.05, value)) %>% 
            left_join(gene_annot, by = "gene_id") %>%
            mutate(plate.prop = paste(high_prop, plate)) %>% 
            as_tibble
  var.arr <- compute.inter.var(ds_deconv)   # This computes coeff of var. 

  # change directory
  #   setwd(args[1])

  # Correlation # ======================================== #
  pdf("outcomes/figures/corr-coefplot2.pdf", width = 9, height = 9, onefile = T)
  par(mfrow = c(2, 2))
  for (plate in c("DPK", "LIT")) {
    for (k in c(0, 1)) {
      rids <- ds_deconv@rdesc$id[ds_deconv@rdesc$high_prop == k]
      ds <- subset.gct(ds_cor, rid = rids)
      d <- extract.cols(ds, handle = handles, plate = plate)
      cor.coefplot(d, x.scale = c(-0.05, 0.15) + ifelse(plate=="DPK", 0.47, 0.62))
      abline(v = mean(d@mat[,10]), lty = 3)
      title(ylab = "solution's rank (holdout data)")
      title(xlab = "mean genewise correlation with 95% CI")
      title(main = sprintf("%s (%s)",gene.prop.labels[k+1], plate.labels[plate]))
    }
  }
  dev.off()
  system("open outcomes/figures/corr-coefplot2.pdf")
  # ======================================== #

  # AUC # ======================================== #
  pdf("outcomes/figures/auc-bars.pdf", width = 8, height = 4)
  par(las = 1, mfrow = c(1, 1), mar = c(1,2,1,1))
  for(p in c("DPK", "LIT")) {
    d <- speed_acc %>% filter(grepl(p, plate)) %>% arrange(rank)
    plot.auc(d$auc, d$method_short, x.baseline = 0.86)
    title(ylab = "AUC extreme modulations", line = 1)
    title(main = plate.labels[p])
  }
  dev.off()
#   system("open outcomes/figures/auc-bars.pdf")
  # ======================================== #

  # KD # ======================================== #
  pdf("outcomes/figures/kd.pdf", width = 5, height = 5)
  par(mar = c(4,4,1,1))
  plot.kd(speed_acc)
  dev.off()
#   system("open outcomes/figures/kd.pdf")

  # Inter-replicate variation   # ======================================== #
  x.pos <- rep(4, 10)
  ord <- order(apply(var.arr$benchmark,2,mean))
  pdf("outcomes/figures/cv-variation.pdf", width = 8, height = 4, onefile = T)
  for (pert in c("BRD","TRCN")) {
    j <- grep(pert, colnames(var.arr$benchmark)[ord], value = T)
    for (k in c(0,1)) {
      par(las = 1, mfrow = c(1, 1), mar = c(4,4,1,1), ann = F)
      genes <- ds_deconv@rdesc$gene_id[ds_deconv@rdesc$high_prop == k]
      y.scale <- 100*range(sapply(var.arr, function(x) apply(x[, j], 2, mean)))
      x.scale <- c(1, length(j)+5)
      #c(5, 20) 
      plot(NA, NA, ylim = y.scale, xlim = x.scale, axes = F)
      box(bty = "l")
      for(i in c(1:4, 10)) {
        handle <- soln_desc$handle[i]
        x.comp <- 100*apply(var.arr[[handle]][genes, j], 2, mean)
        x.se <- 100*apply(var.arr[[handle]][genes, j], 2, function(x)sqrt(var(x)/length(x)))
        lines(x.comp, col = COL.handle[handle], lwd = 2)
        #segments(x0 = 1:length(x.comp), y0 = x.comp + x.se, y1 = x.comp - x.se, col = i)
        x.label <- paste(1:10, soln_desc$method_short)
        x.label[10] <- "dpeak*"
        text(y = tail(x.comp,1), x=length(x.comp), x.label[i], pos = x.pos[i]
            , col = COL.handle[handle])
        text(y = 0, x = 1:length(j), colnames(var.arr$benchmark)[j], srt = 45)
      }
      title(main = sprintf("%s",gene.prop.labels[k+1]))
      title(ylab = "mean coefficient of variation (%)")
      title(xlab = paste(pert.labels[pert], "perturbagen ID"))
      axis(2)
      axis(1, at = c(1, length(j)))
    }
  }
  dev.off()
#   system("open outcomes/figures/cv-variation.pdf")
  # ======================================== #
  
  # Winner ... 
  pdf("outcomes/figures/cor-winner.pdf", width = 8, height = 8, onefile=T)
  plot.winner(cor_df)
  dev.off()

}

main()
