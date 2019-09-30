# identify the best algorithm per gene
# using the testing data. will then
# combine these algos into an ensemble
# and apply to holdout data

library(data.table)
library(cmapR)

# NB: requires access to Broad FS

# gene-level correlations
corr_ds <- parse.gctx("/cmap/projects/crowdsourcing/contest_03_dpeak/results/UNI_DUO_gene_spearman_correlations_n60x976.gctx")

# fix annots - should update in file
corr_ds@cdesc$algorith[corr_ds@cdesc$algorith=="output"] <- "mvaudel"

# sample meta
annot <- fread("/cmap/projects/crowdsourcing/contest_03_dpeak/results/collated_meta.txt")

# deconv data
deconv_ds <- parse.gctx("/cmap/projects/crowdsourcing/contest_03_dpeak/results/DECONV_all_n24772x976.gctx",
                        cid=annot[pert_plate=="DPK.CP003"]$id)

# annotate
deconv_ds <- annotate.gct(deconv_ds, annot, dim="col")

# subset to just the holdout compound plate
deconv_ds_sub <- subset.gct(deconv_ds, cid=which(deconv_ds@cdesc$pert_plate=="DPK.CP003"))

# consider the compound plate from testing data
# and find the algo w/ best corr. per gene
corr_ds_sub <- subset.gct(corr_ds, cid=which(corr_ds@cdesc$plate=="DPK.CP002_MCF7_24H_X1_B42"))
idx <- apply(corr_ds_sub@mat, 1, which.max)
algo2gene <- split(corr_ds_sub@rid, corr_ds_sub@cdesc$algorith[idx])

# make a list of GCT objects, one per algo
algo2ds <- lapply(names(algo2gene), function(x) {
  tmp <- subset.gct(deconv_ds_sub, rid=algo2gene[[x]],
             cid=which(deconv_ds_sub@cdesc$handle==x))
  tmp@cid <- colnames(tmp@mat) <- tmp@cdesc$pert_well
  tmp@cdesc <- data.frame(id=tmp@cid)
  tmp
})
names(algo2ds) <- names(algo2gene)

# wrap up into a matrix
ensemble_ds <- Reduce(function(x, y) {
  merge.gct(x, y, dim="row", matrix_only = T)
}, algo2ds)
ensemble_ds@cdesc <- data.frame(id=ensemble_ds@cid)

# align with ground truth and compute correlation
gt_ds <- subset.gct(deconv_ds_sub, cid=which(deconv_ds_sub@cdesc$handle=="ground-truth"))
gt_ds@cid <- colnames(gt_ds@mat) <- gt_ds@cdesc$pert_well
gt_ds@cdesc <- data.frame(id=gt_ds@cid)
ensemble_ds <- subset.gct(ensemble_ds, rid=gt_ds@rid, cid=gt_ds@cid)
ensemble_corr <- cor(t(ensemble_ds@mat), t(gt_ds@mat), method="spearman")
(ensemble_corr_agg <- median(diag(ensemble_corr)))
# doesn't result in the best overall correlation
# could perhaps be due to only very small differences b/w
# the best algo per gene and/or performance on testing data
# not being entirely indicative of performance on holdout
# we should test this last possibility


# todo - does picking the best algo based on correlation also improve AUC scores?