library(readr)
library(dplyr)

rootdir <- "~/Documents/CATD-ReliabilityAnalysis/"

mdd_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/mdd_v4_icc_schaefersc1000/"), 
                        pattern = "\\.csv$",
                        full.names = TRUE)
mdd_icc_boots <- do.call(rbind.data.frame, lapply(mdd_files, read_csv))
write_csv(mdd_icc_boots, paste0(rootdir, 'data/output/mdd_icc_boots.csv'))

hv_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/hv_v4_icc_schaefersc1000/"), 
                        pattern = "\\.csv$",
                        full.names = TRUE)

hv_icc_boots <- do.call(rbind.data.frame, lapply(hv_files, read_csv, show_col_types = FALSE))

write_csv(hv_icc_boots, paste0(rootdir, 'data/output/hv_icc_boots.csv'))

all_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/all_v4_icc_schaefersc1000/"), 
                       pattern = "\\.csv$",
                       full.names = TRUE)

all_icc_boots <- do.call(rbind.data.frame, lapply(all_files, read_csv, show_col_types = FALSE))

write_csv(all_icc_boots, paste0(rootdir, 'data/output/all_icc_boots.csv'))