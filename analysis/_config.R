# load require packages 
required_packages <- c("dplyr", "rmarkdown", "ggsci", "RColorBrewer", "cmapR", "data.table")
for (package in required_packages)  {
  ok <- suppressMessages(require(package, character.only = T, warn = F, quiet = T))
  if (!ok) {
    install.packages(package)
    library(package)
  }
}

# options -------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE
                    , message = FALSE
                    , warning = FALSE
                    , cache = FALSE
                    , fig.path = "figures/"
                    )

# path to data -------------------------------------------------------
path_data <- "data"
path_speed_acc <- file.path(path_data, "holdout_score_time.txt")
path_soln_desc <- file.path(path_data, "solution_descriptions.txt")
path_annot <- file.path(path_data, "holdout_sample_annotations.txt")
path_gene_annot <- file.path(path_data, "barcode_to_gene_map.txt")
path_ds_cor <- file.path(path_data, "UNI_DUO_gene_spearman_correlations_holdout_n20x976.gctx")
path_ds_deconv <- file.path(path_data, "DECONV_holdout_n8228x976.gctx")


# Colors and labels -------------------------------------------------------
COL           <- ggsci::pal_npg()(9) #brewer.pal(7, "Paired")  
COL.handle    <- c(brewer.pal(9, "Set1"), "#FF7F00")
plate.labels  <- c(DPK = "compound", LIT = "shRNA")
pert.labels   <- c(BRD = "compound", TRCN = "shRNA")
gene.prop.labels <- c("genes in low proportion"
                    , "genes in high proportion")
