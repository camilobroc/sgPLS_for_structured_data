# Author: Camilo Broc (camilo.broc@univ-pau.fr)
# Description: Create synthetic data for the article "Penalized Partial Least Square applied to structured data"

# Library statements

library(MASS)

# Sources ---------------------------

# source(file = "auxiliar.R")
# source(file = "group_index.R")

generate.data.structure <- function(n.obs, l.param, seed, method.generation,seed.sample.true.C,D){
  
  # The function generates synthetic data
  
  ### Args 
  # n.obs : number of observation in the data.
  # l.param : list of parameters for the synthetic datas.
  # seed : seed to reproduce the data.
  # method.generation : equal to "pb5" for first simulation of the article, equal to "pb6" for second simulation of the article
  # seed.sample.true.C : seed for the sampling of the values of the non null elements of C
  # D : The loading vector C.
  
  set.seed(seed)
  
  # Retrieve parameters
  
  p <- l.param$p
  g <- l.param$g
  q <- l.param$q
  n <- n.obs
  s <- l.param$s
  
  # definition of the groups 
  
  list.group.var <- LG.GenerateIsoGroups(g,p)
  list.group.obs <- LG.GenerateIsoGroups(s,n)
  
  
  H <- rep(0,n)
  for (is in 1:s){
    H[LG.GetInterval(list.group.obs,is)] <- rnorm(n = LG.GetLength(list.group.obs,is),
                                                  mean = 0,
                                                  sd = 1)
  }
  
  # Generation of the noise
  
  noise.X <- matrix(0,ncol=p,nrow=n)
  
  for(ig in 1:g){
    for (is in 1:s){
      pg <- LG.GetLength(list.group.var,ig)
      ns <- LG.GetLength(list.group.obs,is)
      noise.X[LG.GetInterval(list.group.obs,is),
        LG.GetInterval(list.group.var,ig)] <- mvrnorm(n = ns,
                                                      mu = rep(0,pg),
                                                      Sigma = l.param$sigma.E.X*(diag(1-l.param$rho,pg) + matrix(l.param$rho,pg,pg)) )
    }
  }
  
  noise.Y <- matrix(0,ncol=q,nrow=n)
  
  for(ig in 1:g){
    for (is in 1:s){
      ns <- LG.GetLength(list.group.obs,is)
      noise.Y[LG.GetInterval(list.group.obs,is),] <- mvrnorm(n = ns,
                                                          mu = rep(0,q),
                                                          Sigma = l.param$sigma.E.Y*(diag(1-l.param$rho,q) + matrix(l.param$rho,q,q)) )
    }
  }
  set.seed(seed.sample.true.C)
  
  # Generation if the loadings depending on the simulation case. pb5 stands for first simulation of the article and pb6 for the second one.
  
  if (method.generation == "pb5"){
    C <- rep(0,p) 
    true.C <- sample(c(rep(1,15),rep(-1,30),rep(1.5,15)))
    C[c(1:15, 21:35, 41:55,61:75)] <- true.C 
    
    # X 
    
    X <- cbind(H)%*%rbind(C)
  } else if (method.generation == "pb6"){
    C <- matrix(0,ncol=p,nrow=3)
    
    true.C <- sample(c(rep(1,15),rep(-1,30),rep(1.5,15)))
    M.C <- matrix(c(1,1,1,1,
                    -1,-1,-1,-1,
                    0.7,0.7,0.7,0.7),ncol = 4,nrow = 3,byrow = TRUE)
    Ind.C <- NULL
    Ind.C[[1]] <- 1:15
    Ind.C[[2]] <- 21:35
    Ind.C[[3]] <- 41:55
    Ind.C[[4]] <- 61:75
    
    Ind.true.C <- NULL
    Ind.true.C[[1]] <- 1:15
    Ind.true.C[[2]] <- 16:30
    Ind.true.C[[3]] <- 31:45
    Ind.true.C[[4]] <- 46:60
    
    for (icl in 1:3){
      for (icc in 1:4){
        C[icl,Ind.C[[icc]]]<-M.C[icl,icc] * rbind(true.C[Ind.true.C[[icc]]])
      }
    }
    X <- matrix(0,ncol = p,nrow = n.obs)
    for (is in 1:s){
      X[LG.GetInterval(list.group.obs,is),] <-  cbind(H[LG.GetInterval(list.group.obs,is)])%*%rbind(C[is,])
    }
    
   
  }
  Y <- cbind(H)%*%rbind(D)
  
  for (ip in 1:p){
    for (is in 1:s){
      ns <- LG.GetLength(list.group.obs,is)
      X[LG.GetInterval(list.group.obs,is),ip] <- X[LG.GetInterval(list.group.obs,is),ip] *l.param$sigma.Batch.X[is] +
        rep(l.param$mu.Batch.X[is],ns)
    }
  }
  
  for (iq in 1:q){
    for (is in 1:s){
      ns <- LG.GetLength(list.group.obs,is)
      Y[LG.GetInterval(list.group.obs,is),iq] <- Y[LG.GetInterval(list.group.obs,is),iq] *l.param$sigma.Batch.Y[is] +
        rep(l.param$mu.Batch.Y[is],ns)
    }
  }
  X <- X + noise.X
  Y <- Y + noise.Y
  
  # print("end")
  return(list(X=X,Y=Y,list.group.obs =list.group.obs, list.group.var = list.group.var))
  
}

translation.shift <- function(X,s,list.group.obs){
  
  # Transforms a matrix into an equivalent matrix following equation (21) of the article.
  
  ### Args 
  # X : Original matrix
  # s : number of observation sets
  # list.group.obs : vector of size s giving the indices of the last element of each observation set
  
  p <- dim(X)[2]
  
  # Translated the matrix blocks
  
  X.translation <- matrix(0,ncol = dim(X)[2]*s,nrow = dim(X)[1])
  for (is in 1:s){
    X.translation[LG.GetInterval(list.group.obs,is), 1:p + (is-1)*p] <- X[LG.GetInterval(list.group.obs,is),] 
  }
  
  # Reorder the columns of the matrix in order to have all columns that will be penalized together being adjecent.
                          
  X.translation.shift <- matrix(0,ncol = dim(X)[2]*s,nrow = dim(X)[1])
  for (is in 1:s){
    X.translation.shift[, ((1:p)-1)*s + is] <- X.translation[,1:p + (is-1)*p] 
  }
  return(X.translation.shift)                        
}

un.translation.shift <- function(X,s,list.group.obs){
  
  # Transforms back an equivalent matrix into the original matrix (equation (21) of the article).
  
  # Reverse the reordering
  p <- dim(X)[2]/s
  X.un.shift<- matrix(0,ncol = p*s,nrow = dim(X)[1])
  for (is in 1:s){
    X.un.shift[,1:p + (is-1)*p] <-  X[, ((1:p)-1)*s + is] 
  }
  
  # Reverse the translation
  X.un.translation <- matrix(0,ncol = p,nrow = dim(X)[1])
  for (is in 1:s){
    X.un.translation[LG.GetInterval(list.group.obs,is),] <-  X.un.shift[LG.GetInterval(list.group.obs,is), 1:p + (is-1)*p] 
  }
  

  return(X.un.translation)                        
}

