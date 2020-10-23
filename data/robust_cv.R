#
# 
# This script creates a file with coefficient of variation array
#
#
library(cmapR)
library(dplyr, warn = F)
path_data <- "../deconv/data"
path_annot      <- file.path(path_data, "holdout_sample_annotations.txt")
path_ds_deconv  <- file.path(path_data, "DECONV_holdout_n8228x976.gctx")
robust_cv <- function (x) {
  100 * IQR(x) / sum(quantile(x, p = c(0.25, 0.75)))
}
get_fun_by_pert <- function(fun, M, pert_id, handle, sig_id, verbose = T, ...) {
  perts <- factor(pert_id)
  handles <- factor(handle)
  out <- array(dim = c(nrow(M), nlevels(perts), nlevels(handles))
              , dimnames = list(rownames(M), levels(perts), levels(handles)))
  for (h in levels(handles)) {
    for (p in levels(perts)) {
      if (verbose) message(h, "/", p)
      sig_id_sel <- sig_id[pert_id == p & handle == h]
      out[, p, h] <- apply(M[, sig_id_sel], 1, fun, ...)
    }
  }
  return(out)
}

# Main -------------------------------------------------------
annot     <- read.table(path_annot, head = T, sep = "\t")
ds_deconv <- parse.gctx(path_ds_deconv)
ds_info <- left_join(annot, ds_deconv@cdesc)
cv <- get_fun_by_pert(robust_cv, ds_deconv@mat
                  , ds_info$pert_id
                  , ds_info$handle
                  , ds_info$id
                  )

# Save -------------------------------------------------------
saveRDS(cv, file = "robust_cv.rds")
  