library(dplyr)
library(readr)
library(ramify)
library(mgc)
library(ggplot2)
library(Hotelling)

rootdir <- '/home/ccc98/Documents/ReliabilityAnalysis/'
all_subs <- read_csv(paste0(rootdir,'references/rest_df_v21_1.csv'), progress = TRUE)
all_connectomes <- read_csv(paste0(rootdir,'data/processed/connectivity_data/icc_dataframe_v21_1.csv'), progress=TRUE)
nboot <- 1000

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

all_subs <- all_subs %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% ungroup()
all_tomes <- all_connectomes %>% filter(grepl("v1|v4", session)) %>% group_by(subject) %>% filter(length(session) == 2) %>% ungroup()
mdd <- all_subs %>% filter(group == "MDD")
hv <- all_subs %>% filter(group == "HV")

mdd_tomes <- merge(mdd, all_tomes, by=c('subject','session'),all.x = TRUE)
hv_tomes <- merge(hv, all_tomes, by=c('subject','session'),all.x = TRUE)


allsubs_v1 <- filter(all_subs, grepl("v1",session))
mdd_v1 <- filter(mdd, grepl("v1",session))
hv_v1 <- filter(hv, grepl("v1",session))

table(allsubs_v1$sex)

table(mdd_v1$sex)

table(hv_v1$sex)


# ICC ---------------------------------------------------------------------

all_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/all_icc_v21_1_boot1000/"), 
                        pattern = "\\.csv$",
                        full.names = TRUE)

all_icc_boots <- lapply(all_files, read_csv)

all_icc_concat <- do.call(rbind.data.frame, all_icc_boots)
all_icc_summ <- all_icc_concat %>% group_by(measure) %>% dplyr::summarise(icc_mean = mean(ICC), icc_sd = sd(ICC), lower_bound = mean(`lower bound`), upper_bound = mean(`upper bound`))



mdd_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/mdd_icc_v21_1_boot1000/"), 
                        pattern = "\\.csv$",
                        full.names = TRUE)

mdd_icc_boots <- lapply(mdd_files, read_csv)

mdd_icc_concat <- do.call(rbind.data.frame, mdd_icc_boots)
mdd_icc_summ <- mdd_icc_concat %>% group_by(measure) %>% summarise(icc_mean = mean(ICC), icc_sd = sd(ICC), lower_bound = mean(`lower bound`), upper_bound = mean(`upper bound`))

hv_files <- list.files(path = paste0(rootdir, "data/output/icc_boots/hv_icc_v21_1_boot1000/"), 
                       pattern = "\\.csv$",
                       full.names = TRUE)

hv_icc_boots <- lapply(hv_files, read_csv)

hv_icc_concat <- do.call(rbind.data.frame, hv_icc_boots)
hv_icc_summ <- hv_icc_concat %>% group_by(measure) %>% summarise(icc_mean = mean(ICC), icc_sd = sd(ICC), lower_bound = mean(`lower bound`), upper_bound = mean(`upper bound`))

mdd_icc_summ$group <- "mdd"
hv_icc_summ$group <- "hv"


write_csv(all_icc_summ,paste0(rootdir,"data/output/bootstrapped1k_all_icc_means_v21_1.csv"))
write_csv(mdd_icc_concat, paste0(rootdir,"data/output/bootstrapped1k_mdd_iccs_v21_1.csv"))
write_csv(hv_icc_concat, paste0(rootdir,"data/output/bootstrapped1k_hv_iccs_v21_1.csv"))
write_csv(mdd_icc_summ, paste0(rootdir,"data/output/bootstrapped1k_mdd_icc_means_v21_1.csv"))
write_csv(hv_icc_summ, paste0(rootdir,"data/output/bootstrapped1k_hv_icc_means_v21_1.csv"))



ht <- hotelling.test(as.matrix(mdd_icc_summ$icc_mean),as.matrix(hv_icc_summ$icc_mean))
# fingerprinting ----------------------------------------------------------

ses1_all <- as.matrix(all_tomes %>% filter(grepl("v1", session)) %>% ungroup() %>% select(-c(1:2)))
ses2_all <- as.matrix(all_tomes %>% filter(grepl("v4", session)) %>% ungroup() %>% select(-c(1:2)))
ses1_mdd <- as.matrix(mdd_tomes %>% filter(grepl("v1", session)) %>% ungroup() %>% select(c("1_0":"121_120")))
ses2_mdd <- as.matrix(mdd_tomes %>% filter(grepl("v4", session)) %>% ungroup() %>% select(c("1_0":"121_120")))
ses1_hv <- as.matrix(hv_tomes %>% filter(grepl("v1", session)) %>% ungroup() %>% select(c("1_0":"121_120")))
ses2_hv <- as.matrix(hv_tomes %>% filter(grepl("v4", session)) %>% ungroup() %>% select(c("1_0":"121_120")))


write.table(ses1_all, paste0(rootdir,'data/processed/ID_code_files/ses1_all.csv'), col.names=FALSE, row.names=FALSE,sep=',')
write.table(ses2_all, paste0(rootdir,'data/processed/ID_code_files/ses2_all.csv'), col.names=FALSE, row.names=FALSE,sep=',')
write.table(ses1_mdd, paste0(rootdir,'data/processed/ID_code_files/ses1_mdd.csv'), col.names=FALSE, row.names=FALSE,sep=',')
write.table(ses2_mdd, paste0(rootdir,'data/processed/ID_code_files/ses2_mdd.csv'), col.names=FALSE, row.names=FALSE,sep=',')
write.table(ses1_hv, paste0(rootdir,'data/processed/ID_code_files/ses1_hv.csv'), col.names=FALSE, row.names=FALSE,sep=',')
write.table(ses2_hv, paste0(rootdir,'data/processed/ID_code_files/ses2_hv.csv'), col.names=FALSE, row.names=FALSE,sep=',')



get_FI <- function(ses1, ses2){
  cor_mat1 <- matrix(0,nrow=dim(ses2)[1],ncol=dim(ses1)[1])
  cor_mat2 <- matrix(0,nrow=dim(ses1)[1],ncol=dim(ses2)[1])
  for (i in 1:dim(ses2)[1]){
    cor_mat1[i,] <- apply(ses1, 1, function(x) {cor.test(x,ses2[i,],method = "s", exact=FALSE)$estimate})
  }
  for (i in 1:dim(ses1)[1]){
    cor_mat2[i,] <- apply(ses2, 1, function(x) {cor.test(x,ses1[i,],method = "s", exact=FALSE)$estimate})
  }
  
  for (i in 1:dim(cor_mat1)[1]) {
    cor_mat1[i,argmax(cor_mat1)[i]] <- 2
    cor_mat2[i,argmax(cor_mat2)[i]] <- 2
  }
  
  ind1 <- cbind(1:dim(cor_mat1)[1], argmax(cor_mat1))
  ind2 <- cbind(1:dim(cor_mat2)[1], argmax(cor_mat2))
  FI1 <- length(ind1[ind1[,1] == ind1[,2],1])/dim(ind1)[1]
  FI2 <- length(ind2[ind2[,1] == ind2[,2],1])/dim(ind2)[1]
  
  return(c(FI1,FI2))
}

FI_all <- get_FI(ses1_all,ses2_all)
FI_mdd <- get_FI(ses1_mdd,ses2_mdd)
FI_hv <- get_FI(ses1_hv,ses2_hv)

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

