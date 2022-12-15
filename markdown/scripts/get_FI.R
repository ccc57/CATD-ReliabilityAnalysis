# -------------------get_FI---------------------------
# Generates fingerprinting index for any two sessions
# Inputs: 
#                   ses1: n-by-n correlation matrix where n is the number of regions
#                   ses2: same thing but for session 2 (order of sessions does not matter)
#                   continuous: boolean, if TRUE, returns vector of length # of subjects with the ratio of within-
#                        subject correlation strength to mean between-subject correlation strength
#                   null: boolean, if TRUE, shuffles matrices to get null value for fingerprinting
#                   
# Outputs:
#                   If continuous is FALSE, returns FI of ses1-ses2 and ses2-ses1 (or null). Take the mean of these for FI.


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