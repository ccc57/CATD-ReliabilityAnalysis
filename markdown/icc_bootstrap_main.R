# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

#Adjust run_name, MDD filter, and sublist depending on the analysis
#Run slurm command file using HPC


rootdir <- '/home/ccc98/Documents/CATD-ReliabilityAnalysis/'
datadir <- paste0(rootdir, 'data/')
run_name <- 'mdd_v4_icc_basc' #CHANGE
group <- 'mdd'
sess <- 'v4' #set to desired timepoint, change all subsequent references

library(readr)
library(dplyr)


conn_datapath <- paste0(datadir,'processed/connectivity_data/connectome_dataframe_basc.csv')
icc_dat <- read_csv(conn_datapath)
icc_dems <- read_csv(paste0(rootdir,'references/rest_df.csv'), progress = TRUE)
# icc_dat_v4 <- icc_dat %>% filter(grepl("v1", session) | grepl("v4", session))
# write_csv(icc_dat_v4,"/home/ccc98/Documents/CATD-ReliabilityAnalysis/data/processed/connectivity_data/data_for_icc_v4_basc.csv")


#Anyone who became an MDD at 1y is labeled MDD
change_status <- icc_dems %>% group_by(subject) %>% filter(length(unique(group)) > 1)
for (i in unique(change_status$subject)){
  icc_dems[icc_dems$subject == i,]$group <- "MDD"
}

# Change MDD filter for timepoint analysis
mdd <- icc_dems %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v4", session)) %>% ungroup() %>% filter(group == "MDD")
# mdd <- icc_dems %>% filter(grepl("v1|v2", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v1", session)) %>% ungroup() %>% filter(group == "MDD")
hv <- icc_dems %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v1", session)) %>% ungroup() %>% filter(group == "HV")
all <- icc_dems %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% filter(grepl("v1", session)) %>% ungroup()


# sublist <- mdd$subject #CHANGE
sublist <- all$subject
# sublist <- hv$subject
nboot <- 1000
set.seed(1)
bootlist <- replicate(nboot, sample(sublist, length(sublist), replace=T))
bootlist[,1] <- sublist
colnames(bootlist) <- sprintf("n%05d", seq(0, 99999))[1:dim(bootlist)[2]]
write_csv(data.frame(bootlist), paste0(rootdir,'references/bootlist_',group,nboot,'.csv'))

icc_path <- paste0(datadir,'processed/connectivity_data/data_for_icc_v4_basc.csv')
bootfile <- paste0(rootdir,'references/bootlist_',group,nboot,'.csv')
fname <- paste0(run_name,nboot)
measure_cols <- colnames(icc_dat %>% select(-c('subject','session')))

x <- 1
sink(paste0(rootdir,'slurm/cmds/',run_name,'.txt'))

# for(bootn in colnames(bootlist)){
#   cat("module load R/4.1.0-foss-2020b; R -s -f ",
#       rootdir,"markdown/scripts/icc_bootrun.R --args ",
#                        rootdir," ",
#                        icc_path," ",
#                        bootfile," ",
#                        bootn, " ",
#                        "0.05 ",
#                        fname, " ",
#                        sess, "\n", sep="")
# }
cat("module load R/4.2.0-foss-2020b;")
for(measure in measure_cols){
  cat("R -s -f ",
      rootdir,"markdown/scripts/icc_bootrun.R --args ",
                       rootdir," ",
                       icc_path," ",
                       bootfile," ",
                       measure, " ",
                       "0.05 ",
                       fname, " ",
                       sess, ";", sep="")
  if((x %% 10) == 0){
    cat("\nmodule load R/4.2.0-foss-2020b;")
  }
  x <- x + 1
}

sink()





# Analysis ----------------------------------------------------------------

