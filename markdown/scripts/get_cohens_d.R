
get_cohens_d <- function(sesA, sesB, measure){
  pooled_sd <- ((sd(sesA[,measure]) * dim(sesA)[1]) + (sd(sesB[,measure]) * dim(sesB)[1]))/(dim(sesA)[1]+dim(sesB)[1])
  d <- (mean(sesA[,measure]) - mean(sesB[,measure])) / pooled_sd
  return(d)
}