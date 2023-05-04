# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

#Takes connectivity dataframe (without subject/session labels) and returns list of subject discriminability values (rank of within-subject distance)
#Called by RS_MDD_Reliability.Rmd
library(data.table)

get_dist_ranklist <- function(tome_data){
  
  distance_matrix <- as.matrix(dist(tome_data, method='euclidian'))
  distance_matrix[upper.tri(distance_matrix, diag=TRUE)] <- NA
  # distcount <- length(distance_matrix[!is.na(distance_matrix)])
  subject_ranklist <- 1:length(within_dist)
  within_dist <- distance_matrix[cbind(seq(1,(dim(distance_matrix)[1]),2) + 1, 
                                 seq(1,(dim(distance_matrix)[1]),2))]
  distance_matrix[cbind(seq(1,(dim(distance_matrix)[1]),2) + 1, 
                        seq(1,(dim(distance_matrix)[1]),2))] <- NA
  distances <- distance_matrix[!is.na(distance_matrix)]
  for (i in 1:length(within_dist)){
    subject_ranklist[i] <- length(distances[distances < within_dist[i]])
  }

  
  return(subject_ranklist)
  
}