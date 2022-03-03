
run_icc <- function(icc_dat,
                    sub_list = unique(icc_dat&subject), 
                    rootdir = 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/',
                    alph = 0.05){

  icc_df <- data.frame()
  sub_list <- data.frame(sub_list)
  psych_dat <- list()
  
  measure_cols <- colnames(icc_dat %>% select(-c('subject','session')))
  pb <- txtProgressBar(min = 0, max = length(measure_cols), style = 3)
  
  psych_files <- list.files(path = paste0(rootdir, "data/processed/connectivity_data/psych_dat_backup"), 
                            pattern = "\\.csv$",
                            full.names = TRUE)
  suppressMessages(psych_dat <- lapply(psych_files, read_csv, show_col_types = FALSE))
  
  
  names(psych_dat) <- gsub("\\.csv$", "", psych_files)
  for(i in 1:length(psych_dat)) {
    names(psych_dat[[i]])[1] <- 'subject'
    write_csv(psych_dat[[i]], paste0(names(psych_dat[i]),'.csv'))
  }
  
  psych_dat = split(psych_dat, f = psych_dat[,1])
  
  for(x in measure_cols){
    i <- which(x == measure_cols)
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

get_boot <- function(dat){
  sublist <- dat$subject
  set.seed(1)
  bootlist <- sample(sublist, dim(sublist)[1], replace=T)
  
}