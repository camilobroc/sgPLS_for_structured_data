# Author: Camilo Broc (camilo.broc@univ-pau.fr)
# Description: File with base functions for indexes of groups. To be sourced
# The functions are called LG for list group


# Functions for indexes of groups

LG.GenerateIsoGroups <- function(g,p){
  # Generate the vector of indices corresponding to g groups/sets of p/g elements 
  
  ng <- p/g
  if (ng-floor(ng)>0){warning("p/g non integer")}
  else {ng <- as.integer(ng)}
  
  x<-(1:g)*ng
  return(x)
}

LG.GetNumberVar <- function(LG){
  
  # number of variables of the dimension (variables or observation)
  
  x <- LG[length(LG)]
  return(x)
}

LG.GetNumberGroup <- function(LG){
  
  # number of groups/sets of the dimension (variables or observation)
  
  x <- length(LG)
  return(x)
}

LG.GetLength <- function(LG,ig){
  
  # Length of the ig-th group/set
  
  if (ig == 1){return(LG[1])}
  else return(LG[ig]-LG[ig-1])
}

LG.GetLengths <- function(LG){
  
  # Vector of the length of each group/set
  
  res <- LG-c(0,LG[-length(LG)])
  return(res)
}

LG.GetInterval <- function(LG,ig){
  
  # Indices corresponding to the ig-th group/set
  
  if (ig == 1){return(1:LG[1])}
  else return((LG[ig-1]+1):LG[ig])
}

LG.GetIntervals <- function(LG){
  
  # Indices corresponding to the length of each group/set
  
  lapply(1:LG.GetNumberGroup(LG),LG.GetInterval,LG = LG)
}



