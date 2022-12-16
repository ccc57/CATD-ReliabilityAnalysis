rootdir <- '/home/ccc98/Documents/ReliabilityAnalysis/'
run_name <- 'mdd_v2_icc_boot' #CHANGE

library(readr)
library(dplyr)

icc_dat <- read_csv(paste0(rootdir,'data/processed/connectivity_data/icc_dataframe_v21_1.csv'))
icc_dems <- read_csv(paste0(rootdir,'references/rest_df_v21_1_light.csv'), progress = TRUE)

# Changing MDD filter to v3 for v3 analysis
# mdd <- icc_dems %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v1", session)) %>% ungroup() %>% filter(group == "MDD")
mdd <- icc_dems %>% filter(grepl("v1|v2", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v1", session)) %>% ungroup() %>% filter(group == "MDD")
hv <- icc_dems %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v1", session)) %>% ungroup() %>% filter(group == "HV")
all <- icc_dems %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v1", session)) %>% ungroup()


sublist <- mdd$subject #CHANGE
nboot <- 1000
set.seed(1)
bootlist <- replicate(nboot, sample(sublist, length(sublist), replace=T))
bootlist[,1] <- sublist
colnames(bootlist) <- sprintf("n%05d", seq(0, 99999))[1:dim(bootlist)[2]]
write_csv(data.frame(bootlist), paste0(rootdir,'references/bootlist_mdd_v2.csv')) #CHANGE



commands <- list()
icc_path <- paste0(rootdir,'data/processed/connectivity_data/icc_dataframe_v21_1.csv')
bootfile <- paste0(rootdir,'references/bootlist_mdd_v2.csv')
fname <- paste0(run_name,nboot)

sink(paste0(rootdir,'slurm/cmds/',run_name,'.txt'))

for(bootn in colnames(bootlist)){
  cat("module load R/4.1.0-foss-2020b; R -s -f ",
      rootdir,"markdown/scripts/icc_bootrun.R --args ",
                       rootdir," ",
                       icc_path," ",
                       bootfile," ",
                       bootn, " ",
                       "0.05 ",
                       fname, "\n", sep="")
}

sink()



# Analysis ----------------------------------------------------------------

