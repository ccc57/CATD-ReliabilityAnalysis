

get_dist_ranklist <- function(tome_data,level){

  if(level == "subject"){
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
  
  else{
    tome_data <- select(all_tomes, -c(1,2))
    pb <- txtProgressBar(min = 0, max = dim(tome_data)[2], style = 3)
    edge_ranklist <- vector()
    for (i in 1:(dim(tome_data)[2])){
      edge_ranklist <- c(edge_ranklist,discr.stat(select(tome_data, -c(i)), all_tomes$subject))
      setTxtProgressBar(pb, i)
      Sys.sleep(.01)
    }
    
    
    
    
    pb <- txtProgressBar(min = 0, max = dim(edge_distance_matrix)[1]-1, style = 3)
    edge_distance_matrix <- as.matrix(dist(t(tome_data), method='euclidian'))
    edge_ranklist <- vector()
    

    edge_distance_matrix[upper.tri(edge_distance_matrix, diag=TRUE)] <- NA
    distcount <- length(edge_distance_matrix[!is.na(edge_distance_matrix)])
    for (i in 1:(dim(edge_distance_matrix)[1]-1)){
      rank <- edge_distance_matrix[i+1,i]
      edge_ranklist <- c(edge_ranklist, 
                            length(edge_distance_matrix[edge_distance_matrix < rank & !is.na(edge_distance_matrix)])/distcount)
      setTxtProgressBar(pb, i)
      Sys.sleep(.01)
    }
    return(edge_ranklist)
  }
}