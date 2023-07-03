library(dplyr)
library(readr)

all_subs <- read_csv(paste0(rootdir, 'references/all_subjects_df.csv'))
all_subs <- all_subs %>% distinct(subject, old_sess, .keep_all = T)

#Filter to subjects with baseline and 1y
all_subs_included <- all_subs %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% 
  filter(length(session) == 2) %>% ungroup()

#Anyone who became an MDD at 1y is labeled MDD
all_subs_included <- all_subs_included %>% 
  group_by(subject) %>% mutate(group = replace(group, "MDD" %in% group, "MDD"))

#Filter to subjects with baseline and 4mo
all_subs_included_v2 <- all_subs %>% filter(!grepl("o|i", session)) %>% 
  filter(grepl("v1|v2", session)) %>% 
  group_by(subject) %>% 
  filter(length(session) == 2) %>% 
  mutate(med_change = rowSums(across(antidepressants:no_meds)))  %>% 
  filter(length(unique(med_change)) == 1) %>% 
  filter(length(unique(group)) == 1) %>% select(-c(med_change)) %>%
  filter(group == "MDD")

mdd <- all_subs_included %>% filter(group == "MDD")
hv <- all_subs_included %>% filter(group == "HV")

schaefersc_connectomes <- read_csv(paste0(rootdir,'data/processed/connectivity_data/connectome_dataframe_schaefersc.csv'), progress=TRUE)
basc_connectomes <- read_csv(paste0(rootdir,'data/processed/connectivity_data/connectome_dataframe_basc.csv'), progress=TRUE)

mdd_sch_tomes <- merge(mdd, schaefersc_connectomes, 
                       by.x=c('subject','old_sess'), 
                       by.y = c("subject","session"),
                       all.x = TRUE)

hv_sch_tomes <- merge(hv, schaefersc_connectomes, 
                       by.x=c('subject','old_sess'), 
                       by.y = c("subject","session"),
                       all.x = TRUE)

all_sch_tomes <- merge(all_subs_included, schaefersc_connectomes, 
                       by.x=c('subject','old_sess'), 
                       by.y = c("subject","session"),
                       all.x = TRUE)

mdd_schaefersc_null_disc_df <- discr.test.one_sample(select(mdd_sch_tomes, -c("subject","session")), mdd_sch_tomes$subject)
# 0.4985157

hv_schaefersc_null_disc_df <- discr.test.one_sample(select(hv_sch_tomes, -c("subject","session")), hv_sch_tomes$subject)
# 0.5019511

all_schaefersc_null_disc_df <- discr.test.one_sample(select(all_sch_tomes, -c("subject","session")), all_sch_tomes$subject)
# 0.4986485