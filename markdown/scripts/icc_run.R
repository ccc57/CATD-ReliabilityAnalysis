
get_psych <- function(rootdir = 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/'){
  psych_files <- list.files(path = paste0(rootdir, "data/processed/connectivity_data/psych_dat"), 
                            pattern = "\\.csv$",
                            full.names = TRUE)
  if(!exists("psych_dat") || length(psych_dat) != length(psych_files)){
    print("Reading files...")
    suppressMessages(psych_dat <- lapply(psych_files, read_csv, show_col_types = FALSE))
    names(psych_dat) <- gsub("\\.csv$", "", psych_files)
  }
  
  for(i in 1:length(psych_dat)){
    psych_dat[[i]]$v1 <- apply(psych_dat[[i]] %>% select(matches("v1")), 1, FUN=function(x){mean(x,na.rm = TRUE)})
    psych_dat[[i]]$v4 <- apply(psych_dat[[i]] %>% select(matches("v4")), 1, FUN=function(x){mean(x,na.rm = TRUE)})
  }
  return(psych_dat)
}


run_icc <- function(icc_dat,
                    sublist = unique(icc_dat$subject), 
                    rootdir = 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/',
                    alph = 0.05,
                    nboot = 1,
                    psych_dat = list()){

  
  measure_cols <- colnames(icc_dat %>% select(-c('subject','session')))
  
  set.seed(1)
  bootlist <- replicate(nboot, sample(sublist, length(sublist), replace=T))
  bootlist[,1] <- sublist
  
  # for(i in 1:length(psych_dat)) {
  #   names(psych_dat[[i]])[1] <- 'subject'
  #   write_csv(psych_dat[[i]], paste0(names(psych_dat[i]),'.csv'))
  # }
  
  psych_files <- list.files(path = paste0(rootdir, "data/processed/connectivity_data/psych_dat"), 
                            pattern = "\\.csv$",
                            full.names = TRUE)
  if(!exists("psych_dat") || length(psych_dat) != length(psych_files)){
    psych_dat <- get_psych(rootdir)
  }
  
  icc_df <- data.frame()
  pb <- txtProgressBar(min = 0, max = length(measure_cols)*nboot, style = 3)
  for(i in 1:length(measure_cols)){
    psych_df <- data.frame()
    for (s in 1:nboot){
      psych_ss <- psych_dat[[i]]
      psych_ss <- psych_ss[match(bootlist[,s],psych_ss$subject),]
      setTxtProgressBar(pb, i*s)
      Sys.sleep(.01)
      suppressMessages(icc <- ICC(psych_ss[,c("v1","v4")], missing=TRUE,alpha = alph/length(measure_cols))$results)
      icc$perm <- s
      psych_df <-rbind(psych_df, icc %>% filter(type == "ICC2"))
    }
    psych_df$measure <- measure_cols[i]
    icc_df <- rbind(icc_df, psych_df)
  }
  return(icc_df)
}
