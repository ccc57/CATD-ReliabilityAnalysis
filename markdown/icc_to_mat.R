library(readr)
library(dplyr)
library(stats)
library(reshape2)
library(ggplot2)

rootdir <- 'C:/Users/chris/OneDrive - Yale University/Projects/ReliabilityAnalysis/'
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

bs_iccs_all <- read_csv(paste0(rootdir,'data/output/bootstrapped_all_icc_means_v21_1.csv'), progress=TRUE)
bs_iccs_all <- aggregate(bs_iccs_all[,1:5], by=list(bs_iccs_all$measure), FUN=mean)
bs_iccs_all[c("measure")] <- list(NULL)
bs_iccs_all <- bs_iccs_all %>% rename(measure = Group.1)
bs_iccs_mdd <- read_csv(paste0(rootdir,'data/output/bootstrapped_mdd_icc_means_v21_1.csv'), progress=TRUE)
bs_iccs_hv <- read_csv(paste0(rootdir,'data/output/bootstrapped_hv_icc_means_v21_1.csv'), progress=TRUE)


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

contrast_mat <- mdd_mat - hv_mat

plot_mtx(all_mat)
plot_mtx(mdd_mat)
plot_mtx(hv_mat)

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



