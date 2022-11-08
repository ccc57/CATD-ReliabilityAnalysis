library(tidyr)
library(dplyr)
library(readr)
library(psych)

run_icc <- function(icc_dat,
                    bootlist = unique(icc_dat$subject),
                    bootn = 1,
                    rootdir = 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/',
                    alph = 0.05,
                    fname = 'icc_mdd'){
  
  
  measure_cols <- colnames(icc_dat %>% select(-c('subject','session')))
  pb <- txtProgressBar(min = 0, max = length(measure_cols), style = 3)
  psych_df <- data.frame()
  for(measure_col in measure_cols){
    psych_dat <- icc_dat[,c("subject","session",measure_col)]
    psych_dat <- pivot_wider(psych_dat, names_from = session, values_from = measure_col)
    psych_dat <- psych_dat[match(bootlist[,c(bootn)],psych_dat$subject),]
    setTxtProgressBar(pb, which(measure_cols == measure_col))
    Sys.sleep(.01)
    suppressMessages(icc <- ICC(psych_dat[,c("v1","v4")], missing=TRUE,alpha = alph/length(measure_cols))$results)
    icc$perm <- bootn
    psych_df <-rbind(psych_df, icc %>% filter(type == "ICC2"))
  }
  write_csv(psych_df, paste0(rootdir,'data/output/icc_boots/',fname,'_boot_',bootn,'.csv'))
}