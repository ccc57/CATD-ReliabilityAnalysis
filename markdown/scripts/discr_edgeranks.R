# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

# Takes NxM connectivity dataframe where N is number of scans and M is number of edges, edge ID (column number), and subject IDs 
# Returns discriminability with an edge removed, referenced by discr_edgeranks_main.R

library(mgc)
library(tidyr)

discr_edgeranks <- function(data, edge_ix, ids){
  return(discr.stat(select(data, -c(edge_ix)), ids))
}