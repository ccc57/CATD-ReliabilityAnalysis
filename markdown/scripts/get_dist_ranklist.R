# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

#Takes connectivity dataframe (without subject/session labels) and returns list of subject discriminability values (rank of within-subject distance)
#Called by RS_MDD_Reliability.Rmd

get_dist_ranklist <- function(tome_data,level){
  
  distance_matrix <- as.matrix(dist(tome_data, method='euclidian'))
  distance_matrix[upper.tri(distance_matrix, diag=TRUE)] <- NA
  distcount <- length(distance_matrix[!is.na(distance_matrix)])
  subject_ranklist <- vector()
  for (i in seq(1,dim(distance_matrix)[1]-1,2)){
    rank <- distance_matrix[i+1,i]
    subject_ranklist <- c(subject_ranklist, 
                          length(distance_matrix[distance_matrix < rank & !is.na(distance_matrix)])/distcount)
  }
  return(subject_ranklist)
  
}