# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

#Takes N subjects x M edges connectivity matrix for each group and edge name, returns cohen's d effect size
#Called from RS_MDD_Reliability.Rmd


get_cohens_d <- function(A, B){
  pooled_sd <- ((sd(A) * (length(A)-1)) + (sd(B) * (length(B)-1)))/(length(A)+length(B))
  d <- (mean(A) - mean(B)) / pooled_sd
  return(d)
}

#Connectomes
ses1_mdd_schaefersc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/schaefersc_v4/ses1_mdd.csv'), col_names = F))
ses2_mdd_schaefersc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/schaefersc_v4/ses2_mdd.csv'), col_names = F))
ses1_hv_schaefersc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/schaefersc_v4/ses1_hv.csv'), col_names = F))
ses2_hv_schaefersc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/schaefersc_v4/ses2_hv.csv'), col_names = F))

#Get effect sizes using the average of both sessions
ses_mean_mdd_sch <- (ses1_mdd_schaefersc+ses2_mdd_schaefersc)/2
ses_mean_hv_sch <- (ses1_hv_schaefersc+ses2_hv_schaefersc)/2

cohens_d_sch <- matrix(sapply(1:ncol(ses_mean_mdd_sch), 
                              function(x){get_cohens_d(ses_mean_mdd_sch[,x], 
                                                       ses_mean_hv_sch[,x])}))

write.table(cohens_d_sch, paste0(rootdir, 'data/output/cohens_d_schaefersc_v4.csv'), row.names = F)


#Repeat for BASC
ses1_mdd_basc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/basc_v4/ses1_mdd.csv'), col_names = F))
ses2_mdd_basc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/basc_v4/ses2_mdd.csv'), col_names = F))
ses1_hv_basc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/basc_v4/ses1_hv.csv'), col_names = F))
ses2_hv_basc <- as.matrix(read_csv(paste0(rootdir, 'ID_scripts/edge_code/data/basc_v4/ses2_hv.csv'), col_names = F))

#Get effect sizes using the average of both sessions
ses_mean_mdd_basc <- (ses1_mdd_basc+ses2_mdd_basc)/2
ses_mean_hv_basc <- (ses1_hv_basc+ses2_hv_basc)/2

cohens_d_basc <- matrix(sapply(1:ncol(ses_mean_mdd_basc), 
                               function(x){get_cohens_d(ses_mean_mdd_basc[,x], 
                                                        ses_mean_hv_basc[,x])}))

write.table(cohens_d_basc, paste0(rootdir, 'data/output/cohens_d_basc_v4.csv'), row.names = F)
