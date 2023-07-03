library(readr)
library(dplyr)

rootdir <- '~/project/CATD-ReliabilityAnalysis/'

all_subs_behavior <- all_subs %>% 
  mutate(psych_meds = rowSums(select(all_subs, antidepressants:ADHD_medication))) %>% 
  group_by(subject) %>% 
  mutate(mean_mfq = mean(s_mfq_tot), 
         mean_shaps = mean(s_shaps_tot),
         mean_ari1w = mean(s_ari1w_tot),
         mean_scared = mean(s_scared_tot),
         mean_motion = mean(motion_per_tr),
         mean_max_fd = mean(max_fd),
         mean_age = mean(age),
         mean_psychmeds = mean(psych_meds),
         mean_othermeds = mean(other),
         mean_tannergrowth = mean(s_tanner_growth),
         mean_tannerpubic = mean(s_tanner_pubic),
         delta_mfq = s_mfq_tot - lag(s_mfq_tot, default = first(s_mfq_tot)),
         delta_shaps = s_shaps_tot - lag(s_shaps_tot, default = first(s_shaps_tot)),
         delta_ari1w = s_ari1w_tot - lag(s_ari1w_tot, default = first(s_ari1w_tot)),
         delta_scared = s_scared_tot - lag(s_scared_tot, default = first(s_scared_tot)),
         delta_psychmeds = psych_meds - lag(psych_meds, default = first(psych_meds)),
         delta_othermeds = abs(other - lag(other, default = first(other))),
         delta_tannergrowth = 
           abs(s_tanner_growth - lag(s_tanner_growth, default = first(s_tanner_growth))),
         delta_tannerpubic = 
           abs(s_tanner_pubic - lag(s_tanner_pubic, default = first(s_tanner_pubic)))) %>% 
  filter(session == "v4")

#Connectomes
ses1_all_schaefersc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/schaefersc_v4/ses1_all.csv'), col_names = F))
ses2_all_schaefersc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/schaefersc_v4/ses2_all.csv'), col_names = F))
ses1_all_basc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/basc_v4/ses1_all.csv'), col_names = F))
ses2_all_basc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/basc_v4/ses2_all.csv'), col_names = F))


# ICC - within subject variance
within_var_mat <- matrix(mapply(function(x,y) sd(c(x,y)),ses1_all_schaefersc, ses2_all_schaefersc), ncol=ncol(ses1_all_schaefersc))
all_subs_behavior$within_var_sch <- apply(within_var_mat, 1, mean)

within_var_mat <- matrix(mapply(function(x,y) sd(c(x,y)),ses1_all_basc, ses2_all_basc), ncol=ncol(ses1_all_basc))
all_subs_behavior$within_var_basc <- apply(within_var_mat, 1, mean)


# Fingerprinting
source(paste0(rootdir, 'markdown/scripts/get_FI.R'))
all_subs_behavior$FI_continuous_sch <- get_FI(ses1_all_schaefersc, ses2_all_schaefersc, continuous = T) 
all_subs_behavior$FI_continuous_basc <- get_FI(ses1_all_basc, ses2_all_basc, continuous = T) 

#Discriminability
disc_sch <- read_csv(paste0(rootdir, 'data/output/discriminability/all/all_v4_notdSchaefersc_subject_discrim.csv'), col_names = F)
disc_basc <- read_csv(paste0(rootdir, 'data/output/discriminability/all/all_v4_notdBasc122_subject_discrim.csv'), col_names = F)
all_subs_behavior$disc_subject_sch <- unlist(disc_sch)
all_subs_behavior$disc_subject_basc <- unlist(disc_basc)

write_delim(all_subs_behavior, paste0(rootdir, 'data/output/all_subs_behavior.csv'), delim = ',')