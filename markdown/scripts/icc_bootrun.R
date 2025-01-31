# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

# -------------------icc_bootrun---------------------------
# Gets ICC values of a dataset of vectorized connectivity matrices
# Inputs: 
#                   rootdir: Root directory
#                   icc_path: path to connectivity dataframe (N scans x M edges)
#                   bootlist_path: path to bootstrap file generated by icc_bootstrap_main.R
#                   bootn: bootstrap file column number
#                   alph: alpha
#                   fname: run name of slurm commmand file
#                   session: name of second session
#                   
# Outputs:
#                   writes ICC file to rootdir/data/output/icc_boots/fname/icc_boot_bootn.csv
# Called by icc_bootstrap_main.R command file output


library(tidyr)
library(dplyr)
library(readr)
library(psych)
# library(svMisc)

args = commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
rootdir = '~/Documents/CATD-ReliabilityAnalysis/'
icc_path = paste0(rootdir,'data/processed/connectivity_data/data_for_icc_v4.csv')
bootlist_path = paste0(rootdir,'references/bootlist_mdd_v4_icc_schaefersc1000.csv')
# bootn = 'n00000'
measure = "1_0"
alph = 0.05
fname = 'test'
session = 'v4'

rootdir = args[1]
icc_path = args[2]
bootlist_path = args[3]
# bootn = args[4]
measure = args[4]
alph = as.double(args[5])
fname = args[6]
session = args[7]

cat(paste0("rootdir = ", args[1],
             "\nicc_path = ", args[2],
             "\nbootlist_path = ", args[3],
             "\nmeasure = ", args[4],
             # "\nbootn = ", args[4],
             "\nalph = ", args[5],
             "\nfname = ", args[6],
             "\nsession = ",args[7]
            ))


icc_dat <- readr::read_csv(icc_path, col_select = c("subject","session",all_of(measure)))
icc_dat$session[grepl("v1",icc_dat$session)] <- "v1"
icc_dat$session[grepl(session,icc_dat$session)] <- session

bootlist <- read_csv(bootlist_path)

# measure_cols <- colnames(icc_dat %>% select(-c('subject','session')))
# pb <- txtProgressBar(min = 0, max = length(measure_cols), style = 3)
# psych_df <- data.frame(matrix(nrow = length(measure_cols)*2, ncol = 9))
# pb <- txtProgressBar(min = 0, max = length(bootlist), style = 3)
psych_df <- data.frame(matrix(nrow = length(bootlist)*2, ncol = 9))
colnames(psych_df) <- c("type",	"ICC",	"F",	"df1",	"df2",	"p",	"lower bound",	"upper bound",	"measure")

# for(measure_col in measure_cols){
#     psych_dat <- icc_dat[,c("subject","session",measure_col)]
#     psych_dat <- pivot_wider(psych_dat, names_from = session, values_from = measure_col)
#     psych_dat <- merge(bootlist[,bootn], psych_dat,by.x = bootn, by.y = "subject", all.x = TRUE)
#     Sys.sleep(.001)
#     setTxtProgressBar(pb, which(measure_cols == measure_col))
#     suppressMessages(icc <- ICC(psych_dat[,c("v1",session)], missing=TRUE,alpha = alph/length(measure_cols))$results)
#     icc$measure <- measure_col
#     # psych_df <-rbind(psych_df, icc %>% filter(type == "ICC2"))
#     psych_df[which(measure_cols == measure_col):(which(measure_cols == measure_col)+1),] <- icc %>% 
#       filter(type == "ICC2" | type == "ICC3")
# }

iter <- 1
for(boot in colnames(bootlist)){
    # psych_dat <- icc_dat[,c("subject","session",measure_col)]    
    # psych_dat <- pivot_wider(psych_dat, names_from = session, values_from = measure_col)
    psych_dat <- pivot_wider(icc_dat, names_from = session, values_from = measure)
    psych_dat <- merge(bootlist[,boot], psych_dat,by.x = boot, by.y = "subject", all.x = TRUE)
    # Sys.sleep(.001)
    # setTxtProgressBar(pb,which(boot == colnames(bootlist)))
    suppressMessages(icc <- ICC(psych_dat[,c("v1",session)], missing=TRUE,alpha = alph/101928)$results)
    icc$measure <- measure
    # psych_df <-rbind(psych_df, icc %>% filter(type == "ICC2"))
    psych_df[iter:(iter+1),] <- icc %>% filter(type == "ICC2" | type == "ICC3")
    iter <- iter+2
}

CI <- as.matrix(rbind((psych_df[1,"ICC"] - quantile(psych_df[psych_df$type == "ICC2", "ICC"] - psych_df[1,"ICC"], c(.1, .9))), 
                   (psych_df[1,"ICC"] - quantile(psych_df[psych_df$type == "ICC3", "ICC"] - psych_df[2,"ICC"], c(.1, .9)))))


mean_df <- psych_df %>% group_by(type, measure) %>% summarise(ICC = mean(ICC), 
                                                              F = mean(F), 
                                                              df1 = mean(df1), 
                                                              df2 = mean(df2), 
                                                              p = mean(p), 
                                                              "lower bound" = mean(`lower bound`), 
                                                              "upper bound" = mean(`upper bound`))
mean_df <- cbind(mean_df,CI)

system2('mkdir',paste0(" -p ",rootdir,'data/output/icc_boots/',fname))
write_csv(mean_df, paste0(rootdir,'data/output/icc_boots/',fname,'/icc_boot1000_',measure,'.csv'))

