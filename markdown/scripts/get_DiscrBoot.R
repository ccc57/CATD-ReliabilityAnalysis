bootstrap <- function(data){
  subject <- sample(seq(1,dim(data)[1],2),replace=FALSE) #Sample without replacement to avoid inflation
  subject <- subject[1:round(.8*length(subject))] #Reduces sample by factor of .8
  return(c(subject,subject+1))
}

get_DiscrBoot <- function(data){
  return(data[bootstrap(data),])
}


