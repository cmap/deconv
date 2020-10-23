#' TSNE plot
#' =========

#+ include = FALSE
source("_config.R")
#COL <- c("#00AFBB", "#E7B800", "#FC4E07")

# read tsne files and other data
tsne_de     <- fread(file.path(path_data, "holdout_tsne_DE.txt"))
tsne_deconv <- fread(file.path(path_data, "holdout_tsne_DECONV.txt"))
soln_desc   <- fread(file.path(path_data, "solution_descriptions.txt"), na.strings="n/a")

# add a row for ground truth
soln_desc <- rbind(soln_desc, data.table(rank="ground-truth",
                                         handle="ground-truth",
                                         method="ground-truth",
                                         method_short="GT"),
                   use.names=T, fill=T)

# extract handle and plate info
tsne_de[, handle := unlist(strsplit(rid, ":"))[1], .(rid)]
tsne_de[, plate := unlist(strsplit(unlist(strsplit(rid, ":"))[2], "_"))[1], .(rid)]
tsne_de[, pert_type := switch(plate,
                              "DPK.CP003" = "compound",
                              "LITMUS.KD019" = "shRNA"),
        .(rid)]
tsne_deconv[, handle := unlist(strsplit(rid, ":"))[1], .(rid)]
tsne_deconv[, plate := unlist(strsplit(unlist(strsplit(rid, ":"))[2], "_"))[1], .(rid)]
tsne_deconv[, pert_type := switch(plate,
                              "DPK.CP003" = "compound",
                              "LITMUS.KD019" = "shRNA"),
        .(rid)]

# merge with solution desc
tsne_de <- merge(tsne_de, soln_desc, by="handle", all.x=T)
tsne_deconv <- merge(tsne_deconv, soln_desc, by="handle", all.x=T)

# set factor levels
levs <- c("DTR", "GMM", "k-means", "CNN", "other", "GT")
tsne_de$method_short <- factor(tsne_de$method_short,
                               levels=levs)
tsne_deconv$method_short <- factor(tsne_deconv$method_short,
                               levels=levs)

#+ deconv, fig.width = 1.5*6, fig.height = 1.5*4, dpi = 300, dev = "png" 
COL.method <-  brewer.pal(8, "Paired")
#par(mfrow = c(2, 2), mar = c(1,1,4,1))
layout(
matrix(c(1, 1, 2, 2, 3, 3
  , 1, 1, 2, 2, 3, 3
  , 4, 5, 6, 7, 8, 9
  , 10, 11, 12, 13, 14, 0), nrow = 4, ncol = 6, byrow = T)
)
par(mar= c(5,1,1,1) + 0.1)

# 
fig.letters <- letters
lettering <- function() {
  title(list(fig.letters[1], cex = 2), adj = 0)
  fig.letters[-1]
}

# panel (a)
plot(TS1 ~ TS2, tsne_deconv, col = COL[factor(pert_type)]
    , pch = 19, bty = "l", cex = .4, axes = F, ann = F)
legend("topleft", col = COL, pch =19, bty = "n"
      , legend = levels(factor(tsne_deconv$pert_type))
      )
fig.letters <- lettering()
title(xlab = "DECONV data", line = 1)

# panel (b)
plot(TS1 ~ TS2, tsne_deconv, col = COL.method[factor(method_short)]
    , pch = 19, bty = "l", cex = .4, axes = F, ann = F)
legend("topleft", col = COL.method, pch =19, bty = "n"
      , legend = levels(factor(tsne_deconv$method_short))
      )
fig.letters <- lettering()
title(xlab = "DECONV data", line = 1)

# panel (c)
plot(TS1 ~ TS2, tsne_de, col = COL.method[factor(method_short)]
    , pch = 19, bty = "l", cex = .4, axes = F, ann = F)
legend("topleft", col = COL.method, pch =19, bty = "n"
      , legend = levels(factor(tsne_de$method_short))
      )
fig.letters <- lettering()
title(xlab = "DE data", line = 1)

# panel (d)
soln_desc <- soln_desc[order(soln_desc$method_short), ]
handles <- soln_desc$handle 
method <- soln_desc$method_short 
par(mar = c(4,1,1,1))
for (h in seq_along(handles)) {
  plot(TS1 ~ TS2, data = tsne_deconv, subset = handle == handles[h]
    , col = COL.method[factor(method_short)]
    , pch = 19, bty = "l", cex = .4, axes = F, ann = F)
#     box(bty = "o")
  title(xlab = paste(method[h], handles[h], sep = "\n"), line = 2)
#  if (h==5) mtext("TS2", side = 2, cex = 0.75)
#  if (h == 10)  mtext("TS1", side = 1, cex = 0.75) 
  if (h == 1) fig.letters <- lettering()
}

 