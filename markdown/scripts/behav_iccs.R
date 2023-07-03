library(psych)
library(dplyr)
library(tidyr)
library(readr)

all_subs_unfiltered <- read_csv(paste0(rootdir,'references/rest_df.csv'), progress = TRUE)
tanner <- read_csv(paste0(rootdir, 'references/tanner.csv'))
dm <- read_csv(paste0(rootdir,"references/rest_paths_no_dates_race.csv"))

all_subs <- merge(all_subs_unfiltered,
                  dm[,c("subject","session","p_demo_screen_background_race")], 
                  by=c("subject","session"),all.x=TRUE)

all_subs$session[grepl("v1",all_subs$session)] <- "v1"
all_subs$session[grepl("v4",all_subs$session)] <- "v4"
all_subs$session[grepl("v2",all_subs$session)] <- "v2"
all_subs$session[grepl("v3",all_subs$session)] <- "v3"

all_subs <- all_subs %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% 
  filter(length(session) > 1) %>% ungroup()

all_subs <- merge(all_subs, tanner, by.x = c("subject","session"), by.y = c("SDAN", "Task_Visit_Type"), all.x = T)

#Anyone who became an MDD at 1y is labeled MDD
change_status <- all_subs %>% group_by(subject) %>% filter(length(unique(group)) > 1)
for (i in unique(change_status$subject)){
  all_subs[all_subs$subject == i,]$group <- "MDD"
}

all_subs_measures <- all_subs %>% select(c("subject","session","motion_per_tr","max_fd","s_mfq_tot","s_shaps_tot","s_ari1w_tot","s_scared_tot"))


get_behav_icc <- function(dat, q){
  dat <- pivot_wider(dat, names_from = session, values_from = all_of(q))
  return(ICC(dat[,c("v1","v4")], missing=TRUE, alpha = .05/20)$results)
}

behav_icc_df <- data.frame(matrix(0,nrow = length(colnames(all_subs_measures)) - 2, ncol = 9))
iter = 1
for(q in colnames(select(all_subs_measures, -c("subject","session")))){
  icc <- get_behav_icc(select(all_subs_measures, c("subject","session",q)), q)
  behav_icc_df[iter:(iter+1),] <- icc %>% filter(type == "ICC2" | type == "ICC3") %>% cbind(q)
  iter <- iter+2
}

colnames(behav_icc_df) <- c("type",	"ICC",	"F", "df1",	"df2",	"p",	"lower bound",	"upper bound",	"measure")


behav.21 <- behav_icc_df %>% filter(type == "ICC2") %>%
  mutate(`ICC(2,1)` = paste0(round(ICC, 2), " (95% CI: ", round(`lower bound`, 2), ", ", round(`upper bound`, 2), ")")) %>%
  select(c("measure","ICC(2,1)"))

behav.31 <- behav_icc_df %>% filter(type == "ICC3") %>%
  mutate(`ICC(3,1)` = paste0(round(ICC, 2), " (95% CI: ", round(`lower bound`, 2), ", ", round(`upper bound`, 2), ")")) %>%
  select(c("measure","ICC(3,1)")) %>% right_join(behav.21, by = "measure")
  
behav_table <- data.frame(Measure = behav.31$measure, `ICC(2,1)` = behav.31$`ICC(2,1)`, `ICC(3,1)` = behav.31$`ICC(3,1)`)
  
rownames(behav_table) <- c("Motion per TR", "Maximum FD", "S-MFQ","S-SHAPS","S-ARI1w", "S-SCARED")
write.table(behav_table, '~/project/CATD-ReliabilityAnalysis/data/output/behavioral_icc_table.txt')
  
  