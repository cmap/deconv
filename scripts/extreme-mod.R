library(cmapR)
library(dplyr, warn = F)
library(kableExtra)

extreme.mod <- function(x) {
  ifelse(abs(x) > 2, "|DE| > 2", "|DE| < 2")
}

# load data 
paths <- list.files("data", pattern = "DE_[^_]*.gct$", full = T)
de_list <- lapply(paths, function(x) parse.gctx(x)@mat)
names(de_list) <- gsub(".gct", "", gsub(".*DE_","", paths))

# Reorder rows and cols
row_ids <- rownames(de_list$gardn999)
col_ids <- colnames(de_list$gardn999)
de_list$benchmark <- de_list$benchmark[row_ids, col_ids]
de_list$UNI <- de_list$UNI[row_ids, col_ids]

# cross tabulate 
tab1 <- addmargins(table(uni = extreme.mod(de_list$UNI), bench = extreme.mod(de_list$benchmark)))
tab2 <- addmargins(table(uni = extreme.mod(de_list$UNI), winner = extreme.mod(de_list$gardn999)))
 
# Save table
kb <- kable(cbind(tab1, tab2)/1000, big.mark = ",", digits = 1, booktabs = T) %>% 
	add_header_above(c(" " = 1, "Benchmark" = 3, "Winning solution" = 3)) %>%
	add_header_above(c("UNI" = 1, "DUO" = 6)) 
cat(kb, file = "outcomes/tab-extreme-mod.html")
system("open outcomes/tab-extreme-mod.html")


