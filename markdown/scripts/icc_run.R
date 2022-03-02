
run_icc <- function(icc_dat,
                    sub_list = unique(icc_dat&subject), 
                    rootdir = 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/',
                    alph = 0.05){

  icc_df <- data.frame()
  sub_list <- data.frame(sub_list)
  measure_cols <- colnames(icc_dat %>% select(-c('subject','session')))
  pb <- txtProgressBar(min = 0, max = length(measure_cols), style = 3)
  
  for(x in measure_cols){
    i <- which(x == measure_cols)
    psych_dat <- read.csv(paste0(rootdir, 'data/processed/connectivity_data/psych_dat/', substring(x,2), '.csv'))
    psych_dat <- merge(psych_dat, sub_list, by.x = 'X', by.y = 'sub_list')
    psych_dat$v1 <- apply(psych_dat %>% select(matches("v1")), 1, FUN=function(x){mean(x,na.rm = TRUE)})
    psych_dat$v4 <- apply(psych_dat %>% select(matches("v4")), 1, FUN=function(x){mean(x,na.rm = TRUE)})
    setTxtProgressBar(pb, i)
    Sys.sleep(.01)
    suppressMessages(icc <- ICC(psych_dat[,c("v1","v4")], missing=TRUE,alpha = alph/length(measure_cols))$results)
    icc$measure <- substring(x,2)
    icc_df <- rbind(icc_df, icc %>% filter(type == "ICC2"))
  }
  return(icc_df)
}