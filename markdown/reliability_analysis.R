library(dplyr)
library(readr)
library(mgc)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(coin)
library(reshape2)
library(corrplot)
library(Hmisc)
library(psych)
library(sjPlot)
library(PerformanceAnalytics)

#Set this to home directory
rootdir <- '~/project/CATD-ReliabilityAnalysis/'

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


# ICC ---------------------------------------------------------------------

all_schaefersc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/all/all_icc_boots_schaefersc.csv'))
hv_schaefersc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/hv/hv_icc_boots_schaefersc.csv'))
mdd_schaefersc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/mdd/mdd_icc_boots_schaefersc.csv'))

all_basc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/all/all_icc_boots_basc.csv'))
hv_basc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/hv/hv_icc_boots_basc.csv'))
mdd_basc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/mdd/mdd_icc_boots_basc.csv'))

mdd_v2_basc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/mdd/mdd_v2_icc_boots_basc.csv'))
mdd_v2_schaefersc_icc_df <- read_csv(paste0(rootdir, 'data/output/icc/mdd/mdd_v2_icc_boots_schaefersc.csv'))


# FI ----------------------------------------------------------------------

mdd_schaefersc_FI_df <- read_csv(paste0(rootdir, 'data/output/FI/mdd/bootstrapped_FI_mdd_schaefersc_v4.csv'),col_names = F)$X1
hv_schaefersc_FI_df <- read_csv(paste0(rootdir, 'data/output/FI/hv/bootstrapped_FI_hv_schaefersc_v4.csv'),col_names = F)$X1
all_schaefersc_FI_df <- read_csv(paste0(rootdir, 'data/output/FI/all/bootstrapped_FI_all_schaefersc_v4.csv'),col_names = F)$X1

mdd_basc_FI_df <- read_csv(paste0(rootdir, 'data/output/FI/mdd/bootstrapped_FI_mdd_basc_v4.csv'),col_names = F)$X1
hv_basc_FI_df <- read_csv(paste0(rootdir, 'data/output/FI/hv/bootstrapped_FI_hv_basc_v4.csv'),col_names = F)$X1
all_basc_FI_df <- read_csv(paste0(rootdir, 'data/output/FI/all/bootstrapped_FI_all_basc_v4.csv'),col_names = F)$X1

FI_all_null <- 0.0115
FI_mdd_null <- 0.0174
FI_hv_null <- 0.0314
# Discriminability --------------------------------------------------------

mdd_schaefersc_disc_df <- as.numeric(read_csv(paste0(rootdir, 'data/output/discriminability/mdd/mdd_v4_notdSchaefersc_bootstrap_discrim.csv'),col_names = F))
hv_schaefersc_disc_df <- as.numeric(read_csv(paste0(rootdir, 'data/output/discriminability/hv/hv_v4_notdSchaefersc_bootstrap_discrim.csv'),col_names = F))
all_schaefersc_disc_df <- as.numeric(read_csv(paste0(rootdir, 'data/output/discriminability/all/all_v4_notdSchaefersc_bootstrap_discrim.csv'),col_names = F))

mdd_basc_disc <- as.numeric(read_csv(paste0(rootdir, 'data/output/discriminability/mdd/mdd_v4_notdBasc122_bootstrap_discrim.csv'),col_names = F))
hv_basc_disc_df <- as.numeric(read_csv(paste0(rootdir, 'data/output/discriminability/hv/hv_v4_notdBasc122_bootstrap_discrim.csv'),col_names = F))
all_basc_disc_df <- as.numeric(read_csv(paste0(rootdir, 'data/output/discriminability/all/all_v4_notdBasc122_bootstrap_discrim.csv'),col_names = F))

#subject level 

mdd_schaefersc_subj_disc_df <- read_csv(paste0(rootdir, 'data/output/discriminability/mdd/mdd_v4_notdSchaefersc_subject_discrim.csv'),col_names = F)$X1
hv_schaefersc_subj_disc_df <- read_csv(paste0(rootdir, 'data/output/discriminability/hv/hv_v4_notdSchaefersc_subject_discrim.csv'),col_names = F)$X1
all_schaefersc_subj_disc_df <- read_csv(paste0(rootdir, 'data/output/discriminability/all/all_v4_notdSchaefersc_subject_discrim.csv'),col_names = F)$X1

mdd_basc_subj_disc_df <- read_csv(paste0(rootdir, 'data/output/discriminability/mdd/mdd_v4_notdBasc122_subject_discrim.csv'),col_names = F)$X1
hv_basc_subj_disc_df <- read_csv(paste0(rootdir, 'data/output/discriminability/hv/hv_v4_notdBasc122_subject_discrim.csv'),col_names = F)$X1
all_basc_subj_disc_df <- read_csv(paste0(rootdir, 'data/output/discriminability/all/all_v4_notdBasc122_subject_discrim.csv'),col_names = F)$X1

#null

# all_schaefersc_null_disc_df <- discr.test.one_sample(select(schaefersc_connectomes, -c("subject","session")), all_tomes$subject)
# 
# mdd_schaefersc_null_disc_df <- discr.test.one_sample(select(schaefersc_connectomes, -c("subject","session")), mdd_tomes$subject)
# 
# hv_schaefersc_null_disc_df <- discr.test.one_sample(select(schaefersc_connectomes, -c("subject","session")), hv_tomes$subject)


# behavioral measures ------------------------------------------------------

#SI Figure 3, Table 3
all_subs_behavior <- read_csv(paste0(rootdir, 'data/output/all_subs_behavior.csv'))
chart.Correlation(all_subs_behavior_df)


all_subs_behavior_df <- all_subs_behavior[c("within_var_sch","within_var_basc",
                                            "FI_continuous_sch","FI_continuous_basc",
                                            "disc_subject_sch","disc_subject_basc",
                                            "mean_mfq", "mean_shaps", "mean_ari1w", 
                                            "mean_scared", "mean_motion", "mean_max_fd", 
                                            "mean_age", "mean_psychmeds","mean_othermeds",
                                            "mean_tannergrowth", "mean_tannerpubic",
                                            "delta_mfq", "delta_shaps", "delta_ari1w", 
                                            "delta_scared", "delta_psychmeds", 
                                            "delta_othermeds", "delta_tannergrowth", 
                                            "delta_tannerpubic")]

attr(all_subs_behavior_df$within_var_sch, "label") <- "Within-Subject Variance (SCH)"
attr(all_subs_behavior_df$within_var_basc, "label") <- "Within-Subject Variance (BASC)"
attr(all_subs_behavior_df$FI_continuous_sch, "label") <- "Fingerprinting (SCH)"
attr(all_subs_behavior_df$FI_continuous_basc, "label") <- "Fingerprinting (BASC)"
attr(all_subs_behavior_df$disc_subject_sch, "label") <- "Participant Discriminability (SCH)"
attr(all_subs_behavior_df$disc_subject_basc, "label") <- "Participant Discriminability (BASC)"
attr(all_subs_behavior_df$mean_mfq, "label") <- "MFQ (mean)"
attr(all_subs_behavior_df$mean_shaps, "label") <- "SHAPS (mean)"
attr(all_subs_behavior_df$mean_ari1w, "label") <- "ARI (mean)"
attr(all_subs_behavior_df$mean_scared, "label") <- "SCARED (mean)"
attr(all_subs_behavior_df$mean_motion, "label") <- "Motion"
attr(all_subs_behavior_df$mean_max_fd, "label") <- "Max Framewise Displacement"
attr(all_subs_behavior_df$mean_age, "label") <- "Age"
attr(all_subs_behavior_df$mean_psychmeds, "label") <- "Meds-Psychoactive (mean)"
attr(all_subs_behavior_df$mean_othermeds, "label") <- "Meds-Other (mean)"
attr(all_subs_behavior_df$delta_mfq, "label") <- "MFQ (change)"
attr(all_subs_behavior_df$delta_shaps, "label") <- "SHAPS (change)"
attr(all_subs_behavior_df$delta_ari1w, "label") <- "ARI (change)"
attr(all_subs_behavior_df$delta_scared, "label") <- "SCARED (change)"
attr(all_subs_behavior_df$delta_psychmeds, "label") <- "Meds-Psychoactive (change)"
attr(all_subs_behavior_df$delta_othermeds, "label") <- "Meds-Other (change)"
attr(all_subs_behavior_df$mean_tannerpubic, "label") <- "Tanner - Pubic (mean)"
attr(all_subs_behavior_df$mean_tannergrowth, "label") <- "Tanner - Growth (mean)"
attr(all_subs_behavior_df$delta_tannerpubic, "label") <- "Tanner - Pubic (change)"
attr(all_subs_behavior_df$delta_tannergrowth, "label") <- "Tanner - Growth (change)"


#SI table 3
tab_corr(all_subs_behavior_df, fade.ns=F,
         corr.method = "spearman",
         triangle = "l",
         use.viewer = TRUE,
         file = paste0(rootdir,"figures/behavior_schaefer_table.doc"))

tab_corr(all_subs_behavior_df, fade.ns=F,
         corr.method = "spearman",
         triangle = "l",
         use.viewer = TRUE,
         file = paste0(rootdir,"figures/behavior_basc_table.doc"))

#SI Figure 3
cormat <- rcorr(as.matrix(all_subs_behavior_df,type="spearman"))
cormat_sch <- cormat$r[c("within_var_sch","FI_continuous_sch","disc_subject_sch"),
                     c("mean_mfq", "mean_shaps", "mean_ari1w", "mean_scared", 
                       "mean_motion", "mean_max_fd", "mean_age",  
                       "mean_psychmeds","mean_othermeds",
                       "mean_tannerpubic","mean_tannergrowth",
                       "delta_mfq", "delta_shaps", 
                       "delta_ari1w", "delta_scared", 
                       "delta_psychmeds", "delta_othermeds",
                       "delta_tannerpubic", "delta_tannergrowth")]
rownames(cormat_sch) <- c("Within-Subject Variance", "Fingerprinting", "Within-Subject Distance")

cormat_sch_p <- cormat$P[c("within_var_sch","FI_continuous_sch","disc_subject_sch"),
                     c("mean_mfq", "mean_shaps", "mean_ari1w", "mean_scared", 
                       "mean_motion", "mean_max_fd", "mean_age",  
                       "mean_psychmeds","mean_othermeds",
                       "mean_tannerpubic","mean_tannergrowth",
                       "delta_mfq", "delta_shaps", 
                       "delta_ari1w", "delta_scared", 
                       "delta_psychmeds", "delta_othermeds",
                       "delta_tannerpubic", "delta_tannergrowth")]
rownames(cormat_sch_p) <- c("Within-Subject Variance", "Fingerprinting", "Within-Subject Distance")

cormat_basc <- cormat$r[c("within_var_basc","FI_continuous_basc","disc_subject_basc"),
                        c("mean_mfq", "mean_shaps", "mean_ari1w", "mean_scared", 
                          "mean_motion", "mean_max_fd", "mean_age",  
                          "mean_psychmeds","mean_othermeds",
                          "mean_tannerpubic","mean_tannergrowth",
                          "delta_mfq", "delta_shaps", 
                          "delta_ari1w", "delta_scared", 
                          "delta_psychmeds", "delta_othermeds",
                          "delta_tannerpubic", "delta_tannergrowth")]
rownames(cormat_basc) <- c("Within-Subject Variance", "Fingerprinting", "Within-Subject Distance")

cormat_basc_p <- cormat$P[c("within_var_basc","FI_continuous_basc","disc_subject_basc"),
                          c("mean_mfq", "mean_shaps", "mean_ari1w", "mean_scared", 
                            "mean_motion", "mean_max_fd", "mean_age",  
                            "mean_psychmeds","mean_othermeds",
                            "mean_tannerpubic","mean_tannergrowth",
                            "delta_mfq", "delta_shaps", 
                            "delta_ari1w", "delta_scared", 
                            "delta_psychmeds", "delta_othermeds",
                            "delta_tannerpubic", "delta_tannergrowth")]
rownames(cormat_basc_p) <- c("Within-Subject Variance", "Fingerprinting", "Within-Subject Distance")


corrplot(cormat_sch,method="color",p.mat = cormat_sch_p, sig.level = .05/length(cormat_sch))
corrplot(cormat_basc,method="color",p.mat = cormat_basc_p, sig.level = .05/length(cormat_sch))



# effect size -------------------------------------------------------------


#from get_cohens_d.R
cohens_d_sch <- read_csv(paste0(rootdir, "data/output/cohens_d_schaefersc_v4.csv"))


#Getting this from HPC - zscored discriminability with each edge removed
edge_discriminability_sch <- unlist(as.vector(read_csv(paste0(rootdir,"data/output/discriminability/all/all_v4_notdSchaefersc_edge_discrim.csv"), col_names = F)))

edge_icc_sch <- all_schaefersc_icc_df %>% filter(type == "ICC2") %>% group_by(measure) %>% summarise(mean_ICC = mean(ICC))

#Differential power and group consistency from Horien et al. code
dp_sch <- as.vector(t(read_csv(paste0(rootdir,'ID_scripts/edge_code/outs/schaefersc_v4/'),col_names=FALSE)))
gc_sch <- as.vector(t(read_csv(paste0(rootdir,'ID_scripts/edge_code/outs/schaefersc_v4/'),col_names=FALSE)))


effect_size_db <- data.frame(`Between-Group Cohens D` = cohens_d_sch$V1, "ICC" = edge_icc$mean_ICC)
colnames(effect_size_db) <- c("Between-Group Cohen's D", "ICC")

effect_size_db$group_consistency <- gc
effect_size_db$differential_power <- dp


#SI table 2
tab_corr(effect_size_db, fade.ns=F,
         corr.method = "pearson",
         triangle = "l",
         use.viewer = TRUE,
         file = paste0(rootdir,"figures/effect_size_pearson_table.doc"))

tab_corr(effect_size_db, fade.ns=F,
         corr.method = "spearman",
         triangle = "l",
         use.viewer = TRUE,
         file = paste0(rootdir,"figures/effect_size_spearman_table.doc"))

#Figure 4
ggplot(effect_size_db, aes(x=`Between-Group Cohen's D`, y=ICC)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="blue", fill="#69b3a2", se=TRUE) + 
  xlab("MDD - HV Cohen's d") + ylab("Edge-level ICC") +
  theme(text = element_text(size=20))

ggplot(effect_size_db, aes(x=cohens_d, y=group_consistency)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="blue", fill="#69b3a2", se=TRUE) + 
  xlab("MDD - HV Cohen's d") + ylab("Group Consistency") +
  theme(text = element_text(size=20))
ggplot(effect_size_db, aes(x=cohens_d, y=differential_power)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="blue", fill="#69b3a2", se=TRUE) + 
  xlab("MDD - HV Cohen's d") + ylab("Differential Power") +
  theme(text = element_text(size=20))
ggplot(effect_size_db, aes(x=cohens_d, y=edge_discriminability)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="blue", fill="#69b3a2", se=TRUE) +
  xlab("MDD - HV Cohen's d") + ylab("Edge Discriminability") +
  theme(text = element_text(size=20))


# plots -------------------------------------------------------------------


all_stats <- data.frame(cbind(rbind(mean(all_schaefersc_icc_df$ICC),
                                    mean(mdd_schaefersc_icc_df$ICC),
                                    mean(hv_schaefersc_icc_df$ICC)),
                              rbind(mean(all_schaefersc_FI_df),
                                    mean(mdd_schaefersc_FI_df),
                                    mean(hv_schaefersc_FI_df)), 
                              rbind(FI_all_null,FI_mdd_null,FI_hv_null),
                              rbind(mean(all_schaefersc_disc_df),
                                    mean(mdd_schaefersc_disc_df),
                                    mean(hv_schaefersc_disc_df)),
                              rbind(mean(all_schaefersc_null_disc_df$null),
                                    mean(mdd_schaefersc_null_disc_df$null),
                                    mean(hv_schaefersc_null_disc_df$null)),                     
                              rbind(quantile(all_schaefersc_FI_df,.05),
                                    quantile(mdd_schaefersc_FI_df, .05), 
                                    quantile(hv_schaefersc_FI_df, .05)),
                              rbind(quantile(all_schaefersc_FI_df,.95),
                                    quantile(mdd_schaefersc_FI_df, .95), 
                                    quantile(hv_schaefersc_FI_df, .95)),
                              rbind(quantile(all_schaefersc_disc_df, .05),
                                    quantile(mdd_schaefersc_disc_df, .05), 
                                    quantile(hv_schaefersc_disc_df, .05)),
                              rbind(quantile(all_schaefersc_disc_df,.95),
                                    quantile(mdd_schaefersc_disc_df, .95), 
                                    quantile(hv_schaefersc_disc_df, .95))),
                        row.names = c("All","MDD","HV"))

colnames(all_stats) <- c('ICC','FI','FI_null','Discr',
                         'Discr_null','FI_lower','FI_upper','D_lower','D_upper')
all_stats$Chance_d <- c("Chance Discriminability","Chance Discriminability","Chance Discriminability")
all_stats$Chance_fi <- c("Chance Fingerprinting","Chance Fingerprinting","Chance Fingerprinting")

#Figure 2A Fingerpinting Bar
ggplot(all_stats,aes(x=rownames(all_stats),y=FI,fill=rownames(all_stats)))+
  geom_bar(stat="identity", width=0.5)  +
  geom_errorbar(aes(ymin=FI_lower, ymax=FI_upper, width=.2),
                position=position_dodge(.9)) + 
  geom_segment(aes(x=c(.75,2.75,1.75),y=FI_null,xend=c(1.25,3.25,2.25),yend = FI_null, linetype=Chance_fi),linewidth=1) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) + 
  scale_linetype_manual("Chance",values=c("Chance Fingerprinting"=2)) +
  guides(fill="none", linetype=guide_legend(title=NULL,keywidth = 3, keyheight = 1)) + 
  ylab("Fingerprinting Index")+ xlab("Condition") + 
  scale_color_manual(values = c("#d45757", "#FECDAA","#6785E3"),aesthetics = "fill") + 
  theme_classic() + theme(text = element_text(size=20)) + 
  theme(legend.justification=c(1,1),legend.position=c(1,1),
        text = element_text(size=20))

#Figure 2B Discriminability Bar
ggplot(all_stats,aes(x=rownames(all_stats),y=Discr,fill=rownames(all_stats)))+ 
  geom_bar(stat="identity", width=0.5)+coord_cartesian(ylim=c(.49,1.05)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
  geom_errorbar(aes(ymin=D_lower, ymax=D_upper, width=.2),
                position=position_dodge(.9)) + 
  geom_segment(aes(x=c(.75,2.75,1.75),
                   y=Discr_null,xend=c(1.25,3.25,2.25),
                   yend = Discr_null,linetype=Chance_d),size=1) + 
  scale_linetype_manual("Chance",values=c("Chance Discriminability"=2)) + 
  guides(fill="none", linetype=guide_legend(title=NULL,keywidth = 3, keyheight = 1)) + 
  ylab("Discriminability")+ xlab("Condition") + 
  scale_color_manual(values = c("#d45757", "#FECDAA","#6785E3"),aesthetics = "fill") + 
  theme_classic() + 
  theme(legend.justification=c(1,1),legend.position=c(1,1),
        text = element_text(size=20))

#SI Figure 1 df


icc_df_for_violin <- pivot_longer(bootstrapped_icc_data,cols=c(1:3),names_to = "Condition",
                                  values_to = "ICC", names_pattern = "(.*)_icc_means")
icc_df_for_violin[icc_df_for_violin == "mdd"] <- "MDD - 1 Year"
icc_df_for_violin[icc_df_for_violin == "hv"] <- "HV - 1 Year"
icc_df_for_violin[icc_df_for_violin == "mdd_v2"] <- "MDD - 4 Months"


ggplot(icc_df_for_violin, aes(x=Condition, y=ICC, fill=Condition)) + 
  geom_violin(trim = FALSE) + 
  geom_boxplot(width=0.1,fill="white",outlier.shape = NA) +
  scale_color_manual(values = c("#FECDAA", "#CACBE4","#6785E3"),aesthetics = "fill") +
  xlab("Condition") + ylab("Mean ICC") +
  theme_bw() + theme(legend.position="none", text = element_text(size=20)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1))
