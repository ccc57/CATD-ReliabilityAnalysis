library(dplyr)
library(readr)
library(ramify)
library(mgc)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(coin)
library(reshape2)
library(corrplot)
library(Hmisc)
library(psych)


rootdir <- 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/'
all_subs_unfiltered <- read_csv(paste0(rootdir,'references/rest_df_v21_1.csv'), progress = TRUE)
all_connectomes <- read_csv(paste0(rootdir,'data/processed/connectivity_data/icc_dataframe_v21_1.csv'), progress=TRUE)
nboot <- 1000

plot_mtx <- function(Dx, main.title="Distance Matrix", xlab.title="Sample Sorted by Source", ylab.title="Sample Sorted by Source") {
  data <- melt(Dx)
  ggplot(data, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(name="dist(x, y)",
                         colours=c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3"),
                         limits=c(min(Dx), max(Dx))) +
    xlab(xlab.title) +
    ylab(ylab.title) +
    theme_bw() +
    ggtitle(main.title)
}



# Demographics ------------------------------------------------------------

all_subs <- all_subs_unfiltered
all_subs$session[grepl("v1",all_subs$session)] <- "v1"
all_subs$session[grepl("v4",all_subs$session)] <- "v4"

change_status <- all_subs %>% group_by(subject) %>% filter(length(unique(group)) > 1)
for (i in unique(change_status$subject)){
  all_subs[all_subs$subject == i,]$group <- "MDD"
}


all_connectomes$session[grepl("v1",all_connectomes$session)] <- "v1"
all_connectomes$session[grepl("v4",all_connectomes$session)] <- "v4"

all_subs <- all_subs %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) > 1) %>% ungroup()
all_tomes <- all_connectomes %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% ungroup()
mdd <- all_subs %>% filter(group == "MDD")
hv <- all_subs %>% filter(group == "HV")

mdd_tomes <- merge(mdd, all_tomes, by=c('subject','session'),all.x = TRUE)
hv_tomes <- merge(hv, all_tomes, by=c('subject','session'),all.x = TRUE)


all_subs_v1 <- filter(all_subs, grepl("v1",session))
mdd_v1 <- filter(mdd, grepl("v1",session))
hv_v1 <- filter(hv, grepl("v1",session))

nsubs <- dim(all_subs_v1)[1]
nsubs_mdd <- dim(mdd_v1)[1]
nsubs_hv <- dim(hv_v1)[1]

table(all_subs_v1$sex)

table(mdd_v1$sex)

table(hv_v1$sex)


# ICC ---------------------------------------------------------------------

mdd_concat_iccs <-read_csv(paste0(rootdir, 'data/output/bootstrapped1k_mdd_iccs_v21_1.csv'))
hv_concat_iccs <-read_csv(paste0(rootdir, 'data/output/bootstrapped1k_hv_iccs_v21_1.csv'))
all_concat_iccs <-read_csv(paste0(rootdir, 'data/output/bootstrapped1k_all_iccs_v21_1.csv'))

mdd_concat_iccs$group <- "MDD"
hv_concat_iccs$group <- "HV"
# all_iccs <- rbind(mdd_concat_iccs, hv_concat_iccs)

all_iccs_means <- all_concat_iccs %>% group_by(measure, group) %>% summarise(ICC = mean(ICC))

all_iccs_means_wide <- pivot_wider(all_iccs_means, names_from = group, values_from = ICC)
all_iccs_means$group <- as.factor(all_iccs_means$group)

all_iccs_by_measure <- pivot_wider(all_iccs, names_from = measure)

independence_test(all_iccs_means$ICC ~ all_iccs_means$group)
icc_summ_means <- all_iccs_means %>% group_by(group) %>% summarise(mean_ICC = mean(ICC))

icc_summ <- read_csv(paste0(rootdir,'data/output/bootstrapped_all_icc_means_v21_1.csv'))
icc_summ_all_means <- icc_summ %>% group_by(measure) %>% summarise(mean_ICC = mean(icc_mean)) %>% pivot_wider(names_from = measure,values_from = mean_ICC)
all_mean <- mean(icc_summ$icc_mean)

wilcox.test(hv_concat_iccs$ICC, mdd_concat_iccs$ICC, alternative = "two.sided")

mdd_concat_iccs$contrasts <- mdd_concat_iccs$ICC - hv_concat_iccs$ICC
icc_contrasts <- mdd_concat_iccs[order(mdd_concat_iccs$contrasts),]
# mdd_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/mdd_icc_v21_boot100/"), 
#                         pattern = "\\.csv$",
#                         full.names = TRUE)
# 
# mdd_icc_boots <- lapply(mdd_files, read_csv)
# 
# mdd_icc_concat <- do.call(rbind.data.frame, mdd_icc_boots)
# mdd_icc_summ <- mdd_icc_concat %>% group_by(measure) %>% summarise(icc_mean = mean(ICC), icc_sd = sd(ICC))
# 
# hv_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/hv_icc_v21_boot100/"), 
#                        pattern = "\\.csv$",
#                        full.names = TRUE)
# 
# hv_icc_boots <- lapply(hv_files, read_csv)
# 
# hv_icc_concat <- do.call(rbind.data.frame, hv_icc_boots)
# hv_icc_summ <- hv_icc_concat %>% group_by(measure) %>% summarise(icc_mean = mean(ICC), icc_sd = sd(ICC))
# 
# mdd_icc_summ$group <- "mdd"
# hv_icc_summ$group <- "hv"
# 
# icc_summ <- rbind.data.frame(hv_icc_summ, mdd_icc_summ)
# icc_summ <- icc_summ %>% group_by(group) %>% summarise(icc_mean_allreg = mean(icc_mean), icc_sd_allreg = mean(icc_sd))
# 


# fingerprinting ----------------------------------------------------------

ses1_all <- as.matrix(all_tomes %>% filter(grepl("v1", session)) %>% ungroup() %>% select(-c(1:2)))
ses2_all <- as.matrix(all_tomes %>% filter(grepl("v4", session)) %>% ungroup() %>% select(-c(1:2)))
ses1_mdd <- as.matrix(mdd_tomes %>% filter(grepl("v1", session)) %>% ungroup() %>% select(c("1_0":"121_120")))
ses2_mdd <- as.matrix(mdd_tomes %>% filter(grepl("v4", session)) %>% ungroup() %>% select(c("1_0":"121_120")))
ses1_hv <- as.matrix(hv_tomes %>% filter(grepl("v1", session)) %>% ungroup() %>% select(c("1_0":"121_120")))
ses2_hv <- as.matrix(hv_tomes %>% filter(grepl("v4", session)) %>% ungroup() %>% select(c("1_0":"121_120")))


# write.table(ses1_all, paste0(rootdir,'ID_scripts/edge_code/data_v21_1/ses1_all.csv'), col.names=FALSE, row.names=FALSE,sep=',')
# write.table(ses2_all, paste0(rootdir,'ID_scripts/edge_code/data_v21_1/ses2_all.csv'), col.names=FALSE, row.names=FALSE,sep=',')
# write.table(ses1_mdd, paste0(rootdir,'ID_scripts/edge_code/data_v21_1/ses1_mdd.csv'), col.names=FALSE, row.names=FALSE,sep=',')
# write.table(ses2_mdd, paste0(rootdir,'ID_scripts/edge_code/data_v21_1/ses2_mdd.csv'), col.names=FALSE, row.names=FALSE,sep=',')
# write.table(ses1_hv, paste0(rootdir,'ID_scripts/edge_code/data_v21_1/ses1_hv.csv'), col.names=FALSE, row.names=FALSE,sep=',')
# write.table(ses2_hv, paste0(rootdir,'ID_scripts/edge_code/data_v21_1/ses2_hv.csv'), col.names=FALSE, row.names=FALSE,sep=',')



get_FI <- function(ses1, ses2, continuous = FALSE, null = FALSE){
  cor_mat1 <- matrix(0,nrow=dim(ses2)[1],ncol=dim(ses1)[1])
  cor_mat2 <- matrix(0,nrow=dim(ses1)[1],ncol=dim(ses2)[1])
  for (i in 1:dim(ses2)[1]){
    cor_mat1[i,] <- apply(ses1, 1, function(x) {cor.test(x,ses2[i,],method = "s", exact=FALSE)$estimate})
  }
  for (i in 1:dim(ses1)[1]){
    cor_mat2[i,] <- apply(ses2, 1, function(x) {cor.test(x,ses1[i,],method = "s", exact=FALSE)$estimate})
  }
  if (continuous == TRUE) {
    within_corrs1 = diag(cor_mat1)
    btwn_corrs1 = cor_mat1
    diag(btwn_corrs1) = 0
    btwn_corrs_avg1 = rowSums(btwn_corrs1) / (dim(btwn_corrs1)[2] - 1)
    
    within_corrs2 = diag(cor_mat2)
    btwn_corrs2 = cor_mat2
    diag(btwn_corrs2) = 0
    btwn_corrs_avg2 = rowSums(btwn_corrs2) / (dim(btwn_corrs2)[2] - 1)
    return(c(btwn_corrs_avg1,btwn_corrs_avg2))
  }
  else {
    ind1 <- cbind(1:dim(cor_mat1)[1], argmax(cor_mat1))
    ind2 <- cbind(1:dim(cor_mat2)[1], argmax(cor_mat2))
    if (null) {
      ind1[,1] <- sample(ind1[,1], length(ind1[,1]))
      ind2[,1] <- sample(ind2[,1], length(ind2[,1]))
      ind1[,2] <- sample(ind1[,2], length(ind1[,2]))
      ind2[,2] <- sample(ind2[,2], length(ind2[,2]))
    }
    FI1 <- length(ind1[ind1[,1] == ind1[,2],1])/dim(ind1)[1]
    FI2 <- length(ind2[ind2[,1] == ind2[,2],1])/dim(ind2)[1]
    return(c(FI1,FI2))
  }
}

FI_all <- get_FI(ses1_all,ses2_all)
FI_mdd <- get_FI(ses1_mdd,ses2_mdd)
FI_hv <- get_FI(ses1_hv,ses2_hv)

FI_all_null <- get_FI(ses1_all,ses2_all, null=TRUE)
FI_mdd_null <- get_FI(ses1_mdd,ses2_mdd, null=TRUE)
FI_hv_null <- get_FI(ses1_hv,ses2_hv, null=TRUE)

poisson.test(35, alternative = "two.sided")

FI_df <- data.frame(rbind(round(nsubs_mdd*mean(FI_mdd)),rbind(round(nsubs_hv*mean(FI_hv)))))
colnames(FI_df) <- "matched"
FI_df$unmatched <- rbind(nsubs_mdd - FI_df$matched[1], nsubs_hv - FI_df$matched[2])
rownames(FI_df) <- c("mdd","hv")

chi <- chisq.test(FI_df,correct=FALSE)

# Discriminability --------------------------------------------------------

corr_dat <- select(all_tomes, c("1_0":"121_120"))
subs <- all_tomes$subject
discr <- discr.test.one_sample(corr_dat, subs)


corr_dat <- select(mdd_tomes, c("1_0":"121_120"))
subs <- mdd_tomes$subject
discr_mdd <- discr.test.one_sample(corr_dat, subs)


corr_dat <- select(hv_tomes, c("1_0":"121_120"))
subs <- hv_tomes$subject
discr_hv <- discr.test.one_sample(corr_dat, subs)



# Plots -------------------------------------------------------------------

all_stats <- data.frame(cbind(rbind(mean(icc_summ_means$mean_ICC),icc_summ_means$mean_ICC[2],icc_summ_means$mean_ICC[1]),
                              rbind(mean(FI_all),mean(FI_mdd),mean(FI_hv)), 
                              rbind(mean(FI_all_null),mean(FI_mdd_null),mean(FI_hv_null)),
                              rbind(discr$stat,discr_mdd$stat,discr_hv$stat),
                              rbind(mean(discr$null),mean(discr_mdd$null),mean(discr_hv$null))),
                        row.names = c('All','MDD','HV'))
colnames(all_stats) <- c('ICC','FI','FI_null','Discr', 'Discr_null')
all_stats$sig_pos <- all_stats$Discr + .1
all_stats$sig_pos_fi <- all_stats$FI + .1


ggplot(all_stats,aes(x=rownames(all_stats),y=FI,fill=rownames(all_stats)))+geom_bar(stat="identity", width=0.5) + ylim(0, 1) + xlab("Condition") +
  ylab("Fingerprinting Index") + scale_color_manual(values = c("#C96468", "#FECDAA","#2A2782"),aesthetics = "fill") +
  geom_segment(aes(x=(c(.75,1.75,2.75)),y=FI_null,xend=c(1.25,2.25,3.25),yend = FI_null), size = 1) +
  geom_text(aes(x=rownames(all_stats),y=sig_pos_fi),size=8,label='***') + 
  theme_light() + theme(legend.position="none", text = element_text(size=20))

ggplot(all_stats,aes(x=rownames(all_stats),y=Discr,fill=rownames(all_stats)))+geom_bar(stat="identity", width=0.5) + ylim(0, 1) + xlab("Condition") +
  ylab("Discriminability") + scale_color_manual(values = c("#C96468", "#FECDAA","#2A2782"),aesthetics = "fill") +
  geom_segment(aes(x=(c(.75,1.75,2.75)),y=Discr_null,xend=c(1.25,2.25,3.25),yend = Discr_null), size = 1) +
  geom_text(aes(x=rownames(all_stats),y=sig_pos),size=8,label='**') + 
  theme_light() + theme(legend.position="none", text = element_text(size=20)) + coord_flip()

# Chisq
chix <- as.table(rbind(c(18, 10), c(17, 24)))
dimnames(chix) <- list(condition = c("MDD", "HV"),
                       FI = c("Correct","Incorrect"))

chi <- chisq.test(chix)

# pwr.chisq.test(w= sqrt())


# Clinical data -----------------------------------------------------------



all_subs_lr <- all_subs %>% 
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
         delta_mfq = s_mfq_tot - lag(s_mfq_tot, default = first(s_mfq_tot)),
         delta_shaps = s_shaps_tot - lag(s_shaps_tot, default = first(s_shaps_tot)),
         delta_ari1w = s_ari1w_tot - lag(s_ari1w_tot, default = first(s_ari1w_tot)),
         delta_scared = s_scared_tot - lag(s_scared_tot, default = first(s_scared_tot)),
         delta_psychmeds = psych_meds - lag(psych_meds, default = first(psych_meds)),
         delta_othermeds = abs(other - lag(other, default = first(other))))

all_subs_lr_v4 <- filter(all_subs_lr, session == "v4")

# ICC
ses_all <- cbind(ses1_all,ses2_all)
all_subs_lr_v4$within_var <- apply(ses_all,1,sd)


# Fingerprinting
FI_continuous_all <- zeros(nrow = dim(ses1_all)[1], ncol = 2)
FI_continuous_all[,1:2] <- get_FI(ses1_all,ses2_all,continuous=TRUE)
all_subs_lr_v4$FI_continuous <- rowMeans(FI_continuous_all)
all_subs_lr_v4$FI_continuous_1 <- FI_continuous_all[,1]
all_subs_lr_v4$FI_continuous_2 <- FI_continuous_all[,2]

# Discriminability
tome_data <- select(all_tomes,-c(1:2))

distance_matrix <- as.matrix(dist(tome_data, method='euclidian'))
distance_matrix[upper.tri(distance_matrix, diag=TRUE)] <- NA


ranklist <- vector()
for (sub in all_subs_v1$subject){
  rank <- distance_matrix[which(all_subs$subject == sub),which(all_subs$subject == sub)][2,1]
  ranklist <- c(ranklist, 
                length(distance_matrix[distance_matrix < rank & !is.na(distance_matrix)])/length(distance_matrix[!is.na(distance_matrix)]))
}

edge_distance_matrix <- as.matrix(dist(t(tome_data), method='euclidian'))
edge_distance_matrix[upper.tri(edge_distance_matrix, diag=TRUE)] <- NA

edge_ranklist <- vector()
for (edge in colnames(corr_dat)){
  rank <- mean(c(edge_distance_matrix[which(colnames(corr_dat) == edge),], edge_distance_matrix[,which(colnames(corr_dat) == edge)]), na.rm=TRUE)
  edge_ranklist <- c(edge_ranklist, rank)
}

all_subs_lr_v4$dist_rank <- ranklist


cormat <- rcorr(as.matrix(select_if(all_subs_lr_v4, is.numeric)),type="spearman")
cormat_r <- cormat$r[c("within_var","FI_continuous","dist_rank"),c("mean_mfq", "mean_shaps", "mean_ari1w", "mean_scared", "mean_motion", "mean_max_fd", "mean_age", 
                                                                          "mean_psychmeds","mean_othermeds", "delta_mfq", "delta_shaps", "delta_ari1w", "delta_scared", "delta_psychmeds", "delta_othermeds")]
cormat_p <- cormat$P[c("within_var","FI_continuous","dist_rank"),c("mean_mfq", "mean_shaps", "mean_ari1w", "mean_scared", "mean_motion", "mean_max_fd", "mean_age", 
                                                                          "mean_psychmeds","mean_othermeds","delta_mfq", "delta_shaps", "delta_ari1w", "delta_scared", "delta_psychmeds", "delta_othermeds")]


plot(all_subs_lr_v4$within_var, all_subs_lr_v4$mean_othermeds)
plot(all_subs_lr_v4$within_var, all_subs_lr_v4$delta_othermeds)
plot(all_subs_lr_v4$within_var, all_subs_lr_v4$delta_othermeds)

plot(all_subs_lr_v4$mean_scared, all_subs_lr_v4$dist_rank)
plot(all_subs_lr_v4$mean_scared, all_subs_lr_v4$dist_rank)
plot(all_subs_lr_v4$mean_scared, all_subs_lr_v4$dist_rank)
plot(all_subs_lr_v4$delta_shaps, all_subs_lr_v4$dist_rank)
plot(all_subs_lr_v4$delta_scared, all_subs_lr_v4$dist_rank)
plot(all_subs_lr_v4$mean_age, all_subs_lr_v4$FI_continuous_1)


corrplot(cormat_r,method="color")
p_plot <- cormat_p
p_plot[p_plot > (.05/45)] <- -1
corrplot(p_plot,method="color")


FI_motion <- ggplot(all_subs_lr_v4, aes(x=FI_continuous, y=mean_motion)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
  xlab("Fingerprinting Ratio") + ylab("Mean Framewise Displacement") +
  theme(text = element_text(size=20))
FI_motion
FI_maxFD <- ggplot(all_subs_lr_v4, aes(x=FI_continuous, y=max_fd)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
  xlab("Fingerprinting Ratio") + ylab("Maximum Framewise Displacement") +
  theme(text = element_text(size=20))
FI_maxFD


model <- lm(FI_continuous_1 ~ sex + mean_mfq + mean_shaps + mean_ari1w + mean_scared + mean_motion + mean_max_fd + mean_age + 
              delta_mfq + delta_shaps + delta_ari1w + delta_scared + delta_psychmeds + delta_othermeds, data = all_subs_lr_v4)

model_summary <- summary(model)




# Effect size analysis ----------------------------------------------------


ses_mean_mdd <- (ses2_mdd+ses1_mdd)/2
ses_mean_hv <- (ses2_hv + ses1_hv)/2
d <- vector()
d_ses2 <- vector()
d_ses1 <- vector()

for (measure in colnames(ses1_mdd)){
  pooled_sd <- ((sd(ses_mean_mdd[,measure]) * nsubs_mdd) + (sd(ses_mean_hv[,measure]) * nsubs_hv))/nsubs
  pooled_sd_ses2 <- ((sd(ses2_mdd[,measure]) * nsubs_mdd) + (sd(ses2_hv[,measure]) * nsubs_hv))/nsubs
  pooled_sd_ses1 <- ((sd(ses1_mdd[,measure]) * nsubs_mdd) + (sd(ses1_hv[,measure]) * nsubs_hv))/nsubs
  d <- c(d,(mean(ses_mean_mdd[,measure]) - mean(ses_mean_hv[,measure])) / pooled_sd)
  d_ses2 <- c(d_ses2,(mean(ses2_mdd[,measure]) - mean(ses2_hv[,measure])) / pooled_sd_ses2)
  d_ses1 <- c(d_ses1,(mean(ses1_mdd[,measure]) - mean(ses1_hv[,measure])) / pooled_sd_ses1)
}

icc_summ_all_means <- icc_summ_all_means[,colnames(ses1_mdd)]
icc <- as.vector(t(icc_summ_all_means))
dp <- as.vector(t(read_csv(paste0(rootdir,'data/processed/connectivity_data/DP_edges_all_valmat.csv'),col_names=FALSE)))
gc <- as.vector(t(read_csv(paste0(rootdir,'data/processed/connectivity_data/mean_gr_con_edges_valmat.csv'),col_names=FALSE)))

effect_size_db <- data.frame(cbind(d, d_ses2, d_ses1))
effect_size_db$cohens_d <- d
effect_size_db$ICC <- icc
effect_size_db$group_consistency <- gc
effect_size_db$differential_power <- dp
effect_size_db$edge_distance_ranks <- edge_ranklist

effect_size_correlations <- rcorr(as.matrix(effect_size_db), type="spearman")
corrplot(effect_size_correlations$r,method="color")
effect_size_correlations$r[effect_size_correlations$P > (.05/12)] <- 0
corrplot(effect_size_correlations$r[c("ICC","group_consistency","differential_power","edge_distance_ranks"),
                                    c("d","d_ses2","d_ses1")])




icc_gg <- ggplot(effect_size_db, aes(x=cohens_d, y=icc)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
  xlab("MDD > HV Cohen's D") + ylab("Mean Edge ICC") +
  theme(text = element_text(size=20))

gc_gg <- ggplot(effect_size_db, aes(x=cohens_d, y=group_consistency)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
  xlab("MDD > HV Cohen's D") + ylab("Group Consistency") +
  theme(text = element_text(size=20))

dp_gg <- ggplot(effect_size_db, aes(x=cohens_d, y=differential_power)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
  xlab("MDD > HV Cohen's D") + ylab("Differential Power") +
  theme(text = element_text(size=20))

rank_gg <- ggplot(effect_size_db, aes(x=cohens_d, y=edge_distance_ranks)) +
  geom_point(alpha = .2) +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
  xlab("MDD > HV Cohen's D") + ylab("Edge Distance Rank") +
  theme(text = element_text(size=20))

icc_gg
gc_gg
dp_gg
rank_gg



# Discr-to-ICC ------------------------------------------------------------

ICCtoD <- function(ICC){
  D <- .5 + (1/pi) * atan(ICC/(sqrt((1-ICC)*(ICC+3))))
  return(D)
}
  
DtoICC
#Wolfram
#MDD ICC from D 0.838483
#HV ICC from D 0.82288

ICCtoD(all_iccs_group_means$mean_ICC[2])
ICCtoD(0.3083967) #mean all ICC
# MDD Discr from ICC 0.5743976
# HV Discr from ICC 0.5484434
# All Discr from ICC 0.5640193

all_iccs_group_means$mean_ICC[2]
