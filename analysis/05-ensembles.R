require(data.table)
source("_config.R")

# inputs
speed_acc <- fread("data/holdout_score_time.txt", na.strings="n/a")
ens_corr <- fread("data/DPK.CP003_PC3_24H_X1_B42_DECONV_correlation_scores.txt")
ens_auc <- fread("data/DPK.CP003_PC3_24H_X1_B42_DE_auc_scores.txt")

# merge
ens <- merge(ens_corr, ens_auc, by=c("algos", "nalgo"))

# compute total runtime for each algo combo
ens[, sec := {
  these_algos <- unlist(strsplit(algos, ":"))
  sum(speed_acc[handle %in% these_algos]$sec)
}, algos]

# stack
speed_acc$nalgo <- 1
setnames(speed_acc, "handle", "algos")
tmp <- rbind(ens[!grepl("benchmark", algos)], 
             speed_acc[plate=="DPK.CP003_PC3_24H_X1_B42", .(algos, cor, auc, sec, nalgo)],
             fill=T, use.names=T)


#+ ensemble_final, fig.width = 4.5*2, fig.height = 4.5
bp.cor <- boxplot(cor ~ nalgo, data = ens, plot = FALSE)
bp.auc <- boxplot(auc ~ nalgo, data = ens, plot = FALSE)
par(mfrow = c(1, 2), mar = c(4,4,1,1))

# left panel
COL <- c(rgb(0,0,1), rgb(1,0,0))
matplot(t(bp.cor$stats), type = "l", lty = c(3,2,1,2,3), col = COL[1], ann = FALSE, axes = FALSE, ylim = c(0.55, 0.6))
axis(1, at = 1:9, 2:10)
axis(2, at = range(bp.cor$stats), format(range(bp.cor$stats), digits = 2))
box(bty = "l")
title(main = "Spearman rank correlation", xlab = "# algorithms", ylab = expression(rho))
legend("bottomright", lty=1:3, c("50th percentile (median)", "25th/75th percentile", "max/min"), bty ="n")

# right panel 
matplot(t(bp.auc$stats), type = "l", col = COL[2], lty = c(3,2,1,2,3), ann = FALSE, axes = FALSE, ylim = c(0.92, 0.95))
axis(1, at = 1:9, 2:10)
axis(2, at = range(bp.auc$stats), format(range(bp.auc$stats), digits = 2))
box(bty = "l")
title(main = "AUC of extreme modulations", xlab = "# algorithms", ylab = "AUC")
legend("bottomright", lty=1:3, c("50th percentile (median)", "25th/75th percentile", "max/min"), bty ="n")
