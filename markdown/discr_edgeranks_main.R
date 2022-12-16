library(tidyr)
library(readr)
library(foreach)
library(doParallel)

rootdir <- '/home/ccc98/Documents/ReliabilityAnalysis/'
all_connectomes <- read_csv(paste0(rootdir,'data/processed/connectivity_data/icc_dataframe_v21_1.csv'), progress=TRUE)

all_connectomes$session[grepl("v1",all_connectomes$session)] <- "v1"
all_connectomes$session[grepl("v4",all_connectomes$session)] <- "v4"
all_tomes <- all_connectomes %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% ungroup()
tome_data <- all_tomes  %>% select(-c(1:2))

source(paste0(rootdir,"/markdown/scripts/discr_edgeranks.R"))

numCores <- detectCores()
numCores
registerDoParallel(numCores)


edge_ranklist <- foreach(i=1:dim(tome_data)[2], .combine=c) %dopar% {
  discr_edgeranks(tome_data,i,all_tomes$subject)$discr
}

write_csv(data.frame(zscore), paste0(rootdir,"data/output/edge_ranklist.csv"))
