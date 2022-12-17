# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'
# Converts ICCs to ICC connectomes


library(readr)
library(dplyr)
library(stats)
library(reshape2)
library(ggplot2)

rootdir <- 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/'
plot_mtx <- function(Dx, xlab.title="Region", ylab.title="Region") {
  data <- melt(Dx)
  ggplot(data, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(name="ICC",
                         colours=(c("white", "#d45757")),
                         limits=c(min(Dx), max(Dx))) +
    xlab(xlab.title) +
    ylab(ylab.title) +
    theme_bw() + theme(aspect.ratio=1) + scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
}

bs_iccs_all <- read_csv(paste0(rootdir,'data/output/bootstrapped1k_all_icc_means_v21_1.csv'), progress=TRUE)
bs_iccs_all <- aggregate(bs_iccs_all[,1:5], by=list(bs_iccs_all$measure), FUN=mean)
bs_iccs_all[c("measure")] <- list(NULL)
bs_iccs_all <- bs_iccs_all %>% rename(measure = Group.1)
bs_iccs_mdd <- read_csv(paste0(rootdir,'data/output/bootstrapped1k_mdd_icc_means_v21_1.csv'), progress=TRUE)
bs_iccs_hv <- read_csv(paste0(rootdir,'data/output/bootstrapped1k_hv_icc_means_v21_1.csv'), progress=TRUE)
mdd_iccs_v2 <- read_csv(paste0(rootdir, 'data/output/mdd_iccs_v2.csv'))
mdd_iccs_v2$icc_mean <- mdd_iccs_v2$ICC
mdd_iccs_v3 <- read_csv(paste0(rootdir, 'data/output/mdd_iccs_v3.csv'))
mdd_iccs_v3$icc_mean <- mdd_iccs_v3$ICC

ROI <- data.frame(do.call(rbind, strsplit(as.vector(bs_iccs_all$measure), split = "_")))
names(ROI) <- c("ROI_1","ROI_2")
ROI$ROI_1 <- as.integer(ROI$ROI_1) + 1
ROI$ROI_2 <- as.integer(ROI$ROI_2) + 1
bs_iccs_all <- cbind(bs_iccs_all, ROI)

diag <- data.frame(cbind(matrix(data=NA,nrow=max(ROI)),1:max(ROI),1:max(ROI)))
colnames(diag) <- c("icc_mean","ROI_1","ROI_2")

bs_iccs_all <- select(bs_iccs_all, c("icc_mean","ROI_1","ROI_2"))
flip <- bs_iccs_all
flip$ROI_1 <- bs_iccs_all$ROI_2
flip$ROI_2 <- bs_iccs_all$ROI_1


bs_iccs_all <- rbind(bs_iccs_all, flip)
bs_iccs_all <- rbind(bs_iccs_all, diag)

bs_iccs_all <- bs_iccs_all[order(bs_iccs_all$ROI_1,bs_iccs_all$ROI_2),]
all_mat <- matrix(data=bs_iccs_all$icc_mean,nrow=max(ROI),ncol=max(ROI))


bs_iccs_mdd <- cbind(bs_iccs_mdd, ROI)
flip$icc_mean <- bs_iccs_mdd$icc_mean
bs_iccs_mdd <- select(bs_iccs_mdd, c("icc_mean","ROI_1","ROI_2"))
bs_iccs_mdd <- rbind(bs_iccs_mdd, flip)
bs_iccs_mdd <- rbind(bs_iccs_mdd, diag)

bs_iccs_mdd <- bs_iccs_mdd[order(bs_iccs_mdd$ROI_1,bs_iccs_mdd$ROI_2),]
mdd_mat <- matrix(data=bs_iccs_mdd$icc_mean,nrow=max(ROI),ncol=max(ROI))

bs_iccs_hv <- cbind(bs_iccs_hv, ROI)
flip$icc_mean <- bs_iccs_hv$icc_mean
bs_iccs_hv <- select(bs_iccs_hv, c("icc_mean","ROI_1","ROI_2"))
bs_iccs_hv <- rbind(bs_iccs_hv, flip)
bs_iccs_hv <- rbind(bs_iccs_hv, diag)

bs_iccs_hv <- bs_iccs_hv[order(bs_iccs_hv$ROI_1,bs_iccs_hv$ROI_2),]
hv_mat <- matrix(data=bs_iccs_hv$icc_mean,nrow=max(ROI),ncol=max(ROI))

mdd_iccs_v2 <- cbind(mdd_iccs_v2, ROI)
flip$icc_mean <- mdd_iccs_v2$icc_mean
mdd_iccs_v2 <- select(mdd_iccs_v2, c("icc_mean","ROI_1","ROI_2"))
mdd_iccs_v2 <- rbind(mdd_iccs_v2, flip)
mdd_iccs_v2 <- rbind(mdd_iccs_v2, diag)

mdd_iccs_v2 <- mdd_iccs_v2[order(mdd_iccs_v2$ROI_1,mdd_iccs_v2$ROI_2),]
v2_mat <- matrix(data=mdd_iccs_v2$icc_mean,nrow=max(ROI),ncol=max(ROI))

mdd_iccs_v3 <- cbind(mdd_iccs_v3, ROI)
flip$icc_mean <- mdd_iccs_v3$icc_mean
mdd_iccs_v3 <- select(mdd_iccs_v3, c("icc_mean","ROI_1","ROI_2"))
mdd_iccs_v3 <- rbind(mdd_iccs_v3, flip)
mdd_iccs_v3 <- rbind(mdd_iccs_v3, diag)

mdd_iccs_v3 <- mdd_iccs_v3[order(mdd_iccs_v3$ROI_1,mdd_iccs_v3$ROI_2),]
v3_mat <- matrix(data=mdd_iccs_v3$icc_mean,nrow=max(ROI),ncol=max(ROI))

plot_mtx(all_mat)
plot_mtx(mdd_mat)
plot_mtx(hv_mat)
plot_mtx(v2_mat)
plot_mtx(v3_mat)

plot_mtx(contrast_mat)

write_quant_matrix <- function(M, quant, fname){
  top <- as.matrix(M > quantile(M,probs=1-(quant/100),na.rm=TRUE))+0
  bottom <- as.matrix(M < quantile(M,probs=quant/100,na.rm=TRUE))+0
  top[is.na(top)] <- 0
  bottom[is.na(bottom)] <- 0
  write.table(top,file=paste0(rootdir,"data/output/",fname,"_top_icc_matrix.csv"),row.names=FALSE,col.names=FALSE,sep = " ")
  write.table(bottom,file=paste0(rootdir,"data/output/",fname,"_bottom_icc_matrix.csv"),row.names=FALSE,col.names=FALSE,sep = " ")
  return(quantile(M, 1-(quant/100), na.rm=TRUE))
}

quantile(mdd_mat, .99, na.rm=TRUE)

write_quant_matrix(hv_mat, 1, "hv_1")



plot_mtx(mdd_mat_top)
plot_mtx(mdd_mat_bottom)

write.table(all_mat,file=paste0(rootdir,"data/output/all_icc_matrix.txt"),row.names=FALSE,col.names=FALSE,sep = " ")


write.table(mdd_mat,file=paste0(rootdir,"data/output/mdd_icc_matrix.csv"),row.names=FALSE,col.names=FALSE,sep = " ")
write.table(mdd_mat_top,file=paste0(rootdir,"data/output/mdd_icc_top_1.csv"),row.names=FALSE,col.names=FALSE,sep = " ")
write.table(mdd_mat_bottom,file=paste0(rootdir,"data/output/mdd_icc_bottom_1.csv"),row.names=FALSE,col.names=FALSE,sep = " ")



