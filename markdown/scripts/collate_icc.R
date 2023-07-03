library(readr)
library(dplyr)

rootdir <- "~/project/CATD-ReliabilityAnalysis/"

mdd_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/mdd_v2_icc_schaefersc1000/"), 
                        pattern = "\\.csv$",
                        full.names = TRUE)
mdd_icc_boots <- do.call(rbind.data.frame, lapply(mdd_files, read_csv, show_col_types = FALSE))
write_csv(mdd_icc_boots, paste0(rootdir, 'data/output/icc/mdd/mdd_v2_icc_boots_schaefersc.csv'))

hv_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/hv_v2_icc_basc1000/"), 
                        pattern = "\\.csv$",
                        full.names = TRUE)

hv_icc_boots <- do.call(rbind.data.frame, lapply(hv_files, read_csv, show_col_types = FALSE))

write_csv(hv_icc_boots, paste0(rootdir, 'data/output/icc/hv/hv_v2_icc_boots_basc.csv'))

all_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/all_v2_icc_basc1000/"), 
                       pattern = "\\.csv$",
                       full.names = TRUE)

all_icc_boots <- do.call(rbind.data.frame, lapply(all_files, read_csv, show_col_types = FALSE))

write_csv(all_icc_boots, paste0(rootdir, 'data/output/icc/all/all_v2_icc_boots_basc.csv'))