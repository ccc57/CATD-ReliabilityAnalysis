# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

#Takes N subjects x M edges connectivity matrix for each group and edge name, returns cohen's d effect size
#Called from RS_MDD_Reliability.Rmd


get_cohens_d <- function(sesA, sesB, measure){
  pooled_sd <- ((sd(sesA[,measure]) * dim(sesA)[1]) + (sd(sesB[,measure]) * dim(sesB)[1]))/(dim(sesA)[1]+dim(sesB)[1])
  d <- (mean(sesA[,measure]) - mean(sesB[,measure])) / pooled_sd
  return(d)
}