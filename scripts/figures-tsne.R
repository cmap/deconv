library(data.table)
library(ggplot2)

# read tsne files and other data
tsne_de <- fread("data/holdout_tsne_DE.txt")
tsne_deconv <- fread("data/holdout_tsne_DECONV.txt")
soln_desc <- fread("data/solution_descriptions.txt", na.strings="n/a")

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

# First panel 
p_deconv_color_type <- ggplot(tsne_deconv, aes(x=TS1, y=TS2, color=pert_type)) +
  geom_point(alpha=0.5, size=1) +
  guides(color=guide_legend(override.aes = list(alpha=1, size=1))) +
  scale_color_discrete(name="") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  ggtitle("DECONV colored by perturbagen type")


# Second panel 
p_deconv_color_method <- ggplot(tsne_deconv, aes(x=TS1, y=TS2, color=method_short)) +
  geom_point(alpha=0.5, size=0.5) +
  guides(color=guide_legend(override.aes = list(alpha=1, size=1))) +
  scale_color_discrete(name="") +
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  ggtitle("DECONV colored by algorithm type")

# Third panel 
p_deconv_facet_handle <- p_deconv_color_method + facet_wrap(~method_short+handle)


# Fourth 
p_de <- ggplot(tsne_de, aes(x=TS1, y=TS2, color=method_short)) +
  geom_point(alpha=0.5, size=0.5) +
  guides(color=guide_legend(override.aes = list(alpha=1, size=1))) +
  scale_color_discrete(name="") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  facet_wrap(~method_short+handle) +
  ggtitle("DE colored by algorithm type")



pdf("outcomes/figures/tsne.pdf", height = 7, width = 12)
cowplot::plot_grid(p_deconv_color_type,
                   p_deconv_color_method,
                   p_deconv_facet_handle,
                   p_de,
                   align="hv", ncol=2, nrow=2, rel_heights = c(1, 1.6),
                   labels=letters[1:4])
dev.off()