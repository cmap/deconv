#
source("_config.R")

# Load data sets -------------------------------------------------------
speed_acc   <- read.table(path_speed_acc, na.strings="n/a", head = T, sep = "\t")
soln_desc   <- read.table(path_soln_desc, na.strings="n/a", head = T, sep = "\t")
annot       <- read.table(path_annot, head = T, na.strings = "-666", sep = "\t")
gene_annot  <- read.table(path_gene_annot, colClasses = c("gene_id"="character"), head = T)

# process -------------------------------------------------------
setnames(soln_desc, "method_short", "category")

# table 1 -------------------------------------------------------
knitr::kable(soln_desc) %>% 
cat(file = "tab1.md", fill = T)

# table 2 -------------------------------------------------------
annot %>% 
  filter(grepl("BRD", pert_id), handle == "ground-truth") %>%
  count(pert_iname, pert_id, pert_itime, pert_idose) %>%
  rename(num_replicates = n) %>%
  knitr::kable() %>%
  cat(file = "tab2.md", fill = T)


# table 3 -------------------------------------------------------
annot %>% 
  filter(grepl("TRCN", pert_id), handle == "ground-truth") %>%
  count(pert_iname, pert_id) %>%
  rename(num_replicates = n) %>%
  knitr::kable() %>% 
  cat(file = "tab3.md", fill = T)
  