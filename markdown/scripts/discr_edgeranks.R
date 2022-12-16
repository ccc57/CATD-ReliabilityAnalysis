library(mgc)
library(tidyr)

discr_edgeranks <- function(data, edge_ix, ids){
  return(discr.stat(select(data, -c(edge_ix)), ids))
}