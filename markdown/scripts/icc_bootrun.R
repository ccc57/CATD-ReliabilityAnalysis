
library(tidyr)
library(dplyr)
library(readr)
library(psych)
library(svMisc)

args = commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
rootdir = '~/Documents/ReliabilityAnalysis/'
icc_path = paste0(rootdir,'data/processed/connectivity_data/icc_dataframe_v21_1.csv')
bootlist_path = paste0(rootdir,'references/bootlist_mdd_v21_1.csv')
bootn = 'n00000'
alph = 0.05
fname = 'mdd_icc_v21_1_sublevel'
level = "session"

rootdir = args[1]
icc_path = args[2]
bootlist_path = args[3]
bootn = args[4]
alph = as.double(args[5])
fname = args[6]
# level = args[7]

cat(paste0("rootdir = ", args[1],
             "\nicc_path = ", args[2],
             "\nbootlist_path = ", args[3],
             "\nbootn = ", args[4],
             "\nalph = ", args[5],
             "\nfname = ", args[6],
             "\nlevel = ", args[7]))

icc_dat <- read_csv(icc_path)
icc_dat$session[grepl("v1",icc_dat$session)] <- "v1"
icc_dat$session[grepl("v2",icc_dat$session)] <- "v2" #CHANGE FOR SESSION

bootlist <- read_csv(bootlist_path)

measure_cols <- colnames(icc_dat %>% select(-c('subject','session')))
pb <- txtProgressBar(min = 0, max = length(measure_cols), style = 3)
psych_df <- data.frame()

for(measure_col in measure_cols){
    psych_dat <- icc_dat[,c("subject","session",measure_col)]
    psych_dat <- pivot_wider(psych_dat, names_from = session, values_from = measure_col)
    psych_dat <- merge(bootlist[,bootn], psych_dat,by.x = bootn, by.y = "subject", all.x = TRUE)
    setTxtProgressBar(pb, which(measure_cols == measure_col))
    Sys.sleep(.01)
    suppressMessages(icc <- ICC(psych_dat[,c("v1","v2")], missing=TRUE,alpha = alph/length(measure_cols))$results) 
    icc$measure <- measure_col
    psych_df <-rbind(psych_df, icc %>% filter(type == "ICC2"))
}

system2('mkdir',paste0(" -p ",rootdir,'data/output/icc_boots/',fname))
write_csv(psych_df, paste0(rootdir,'data/output/icc_boots/',fname,'/icc_boot_',bootn,'.csv'))

