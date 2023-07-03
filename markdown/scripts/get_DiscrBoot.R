# title: "Test-retest reliability of functional connectivity in depressed adolescents"
# author: "Chris C. Camp, Stephanie Noble, Dustin Scheinost, Argyris Stringaris, and Dylan M. Nielson"
# date: '2022-12-17'

#Helper function for bootstrapped discriminability
#Takes connectivity dataframe and returns dataframe sampled without replacement at .8 sample size
#Called by RS_MDD_Reliability.Rmd

bootstrap <- function(data){
  subject <- sample(seq(1,dim(data)[1],2),replace=FALSE) #Sample without replacement to avoid inflation
  subject <- subject[1:round(.8*length(subject))] #Reduces sample by factor of .8
  return(c(subject,subject+1))
}

get_DiscrBoot <- function(data){
  return(data[bootstrap(data),])
}


