# Author: Camilo Broc (camilo.broc@univ-pau.fr)
# Description: Application of sparse group PLS on structured data


library(ggplot2)
library(lattice)
library(MASS)
library(mixOmics)
library(Rcpp)

# Sources ---------------------------

# Those functions come from the package "sgspls" from Matthew Sutton 
source(file = "internal.R")
source(file = "internalPLS.R")
source(file = "sgspls.R")
source(file = "predict.R")
source(file = "coef.R")
source(file = "perf.R")
source(file = "tuning.R")
sourceCpp("sgsPLS.cpp")


# Sourcing other files
source(file = "group_index.R")
source(file = "generation.R")
source(file = "sgPLS_structure.R")

{
  # Setting of the parameters.  
  
  # mu.Batch.X : corresponds to mu^(X,B) in the article 
  # sigma.Batch.X : corresponds to lambda^(X,B) in the article
  # mu.Batch.Y :  corresponds to mu^(Y,B) in the article
  # sigma.Batch.Y :  corresponds to lambda^(Y,B) in the article
  # D :  corresponds to D in the article
  # sigma.E.X : corresponds to lambda^(X,E) in the article. This is the noise level.
  # sigma.E.Y :  corresponds to lambda^(Y,E) in the article. This is the noise level.
  
l.param <- list(mu.Batch.X = 0*c(1 ,-0.5 ,-0.5),
                sigma.Batch.X = 1*c(1,0.8,1.5),
                mu.Batch.Y =0*c(1, 0 ,-1),
                sigma.Batch.Y = 1*c(0.6,1.4,1),
                D=c(1,-1,1.5),
                sigma.E.X =2,
                sigma.E.Y =2,
                p = 1000,
                g = 50,
                q = 3,
                s = 3,rho=0.05
                )
  
  # setting for the loop of repetitions
  
  nrep=50
  seed.0 =20
  result <- matrix(0,ncol=6,nrow=4)
  colnames(result)<- c("MSEP:Y1","MSEP:Y2","MSEP:Y3","MSEP:Total","TPR","TD")
  rownames(result)<- c("spls","MINT","sPLSs","sgPLSs")
  
  # Loop of 50 repetitions
  
  for (irep in 1:nrep){
    seed = (irep+ seed.0)
    print(seed)
res.spls <- my.perf.spls(keepX = 60,n.train = 900,n.test = 300, seed = seed, l.param = l.param,method="spls",method.generation = "pb6")
res.mint <- my.perf.spls(keepX = 60,n.train = 900,n.test = 300, seed = seed, l.param = l.param,method="MINT",method.generation = "pb6")
# truc <- system.time(
# res.sgpls <- my.perf.spls(keepX = 4,n.train = 900,n.test = 300, seed = 2, l.param = l.param,method="sgpls",method.generation = "pb5")
# )
res2.spls <- my.perf.spls(keepX = 60,n.train = 900,n.test = 300, seed = seed, l.param = l.param,method="spls_structure",method.generation = "pb6")
truc <- system.time(
  res2.sgpls <- my.perf.spls(keepX = 4,n.train = 900,n.test = 300, seed = seed, l.param = l.param,method="sgpls_structure",method.generation = "pb6")
)


result <- result + (rbind(
  rbind(c(res.spls$sum.square[1:4],res.spls$TPR,res.spls$TD)),
  rbind(c(res.mint$sum.square[1:4],res.mint$TPR,res.mint$TD))
  ,rbind(c(res2.spls$sum.square[1:4],res2.spls$TPR,res2.spls$TD))
  ,rbind(c(res2.sgpls$sum.square[1:4],res2.sgpls$TPR,res2.sgpls$TD))
))
}
result <- result/nrep
result
}
save(result,file="result_second_case.RData")

