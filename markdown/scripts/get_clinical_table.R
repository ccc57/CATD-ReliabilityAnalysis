library(readr)
library(dplyr)

#Set this to home directory
rootdir <- '~/project/CATD-ReliabilityAnalysis/'

all_subs <- read_csv(paste0(rootdir, 'references/all_subjects_df.csv'))
all_subs <- all_subs %>% distinct(subject, old_sess, .keep_all = T)

dm <- read_csv(paste0(rootdir,"references/rest_paths_no_dates_race.csv"))

all_subs_dm <- merge(all_subs, dm[,c("subject","session","p_demo_screen_background_hispanic")],
                     by.x = c("subject","old_sess"), by.y = c("subject","session"),all.x = T)


#Filter to subjects with baseline and 1y
all_subs_included_v1 <- all_subs %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% 
  filter(length(session) == 2) %>% ungroup() %>% mutate(group_time = "v4")

#Filter to subjects with baseline and 4mo
all_subs_included_v2 <- all_subs %>% filter(!grepl("o|i", session)) %>% 
  filter(grepl("v1|v2", session)) %>% 
  group_by(subject) %>% 
  filter(length(session) == 2) %>% 
  mutate(med_change = rowSums(across(antidepressants:no_meds)))  %>% 
  filter(length(unique(med_change)) == 1) %>% 
  filter(length(unique(group)) == 1) %>% select(-c(med_change)) %>% 
  filter(group == "MDD") %>%  mutate(group_time = "v2")

all_subs_included <- rbind(all_subs_included_v1, all_subs_included_v2)

all_subs_excluded <- all_subs %>% filter(!(subject %in% all_subs_included$subject)) %>% 
  filter(session == "v1") %>% mutate(group_time = "excluded")

all_subs_table <- rbind(all_subs_included, all_subs_excluded)

#Anyone who became an MDD at 1y is labeled MDD
all_subs_table[all_subs_table$group_time != "v2",] <- all_subs_table[all_subs_table$group_time != "v2",] %>% 
  group_by(subject) %>% mutate(group = replace(group, "MDD" %in% group, "MDD"))

tbl <- all_subs_table %>% group_by(group_time, group, session) %>% 
  summarise(n = length(unique(subject)), 
            Males = sum(sex == "M"),
            Females = sum(sex == "F"),
            White = sum(p_demo_screen_background_race == "white",na.rm = T),
            Black = sum(p_demo_screen_background_race == "black", na.rm = T),
            Asian = sum(p_demo_screen_background_race == "asian", na.rm = T),
            Pacific_Islander = sum(p_demo_screen_background_race == "pacific", na.rm = T),
            Biracial = sum(p_demo_screen_background_race == "biracial", na.rm = T),
            No_race = sum(is.na(p_demo_screen_background_race) | p_demo_screen_background_race == "decline_answer"),
            Age = mean(age)/12,
            MFQ = mean(s_mfq_tot, na.rm = T),
            MFQ_c = sum(!is.na(s_mfq_tot)),
            SHAPS = mean(s_shaps_tot, na.rm = T),
            SHAPS_c = sum(!is.na(s_shaps_tot)),
            ARI = mean(s_ari1w_tot, na.rm = T),
            ARI_c = sum(!is.na(s_ari1w_tot)),
            SCARED = mean(s_scared_tot, na.rm = T),
            SCARED_c = sum(!is.na(s_scared_tot)))

write.table(tbl, paste0(rootdir, 'figures/demo_table.csv'), row.names = F, sep = ',')
