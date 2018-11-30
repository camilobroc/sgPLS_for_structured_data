# Author: Camilo Broc (camilo.broc@univ-pau.fr)
# Description: Performs sPLS, MINT, sgPLS, sPLS for structured data and sgPLS for structured data and evaluate their performance.

standardization <- function(x){return((x-mean(x))/sd(x))}

standardize.by.study <- function(X,Y,l.obs){
  
  # Standardize the study acccording to the observation sets
  
  ### Args
  # X,Y : the two studied matrices 
  # l.obs : vector corresponding to the observation sets.
  
  ### Out
  # X,Y : standardized matrices
  # mean.X, sd.X : (s x dim(X)[1]) matrices of resp. mean and standard deviation by which X is standardized. each row correspond to a different observation set
  # mean.Y, sd.Y : same for Y
  
  s <- length(l.obs)

  X.std.p.s <- X
  mean.X <- NULL
  sd.X <- NULL
  for (is in 1:s){
    mean.X <- rbind(mean.X, apply(X=X[LG.GetInterval(l.obs,is),],MARGIN = 2,FUN=mean))
    sd.X <- rbind(sd.X, apply(X=X[LG.GetInterval(l.obs,is),],MARGIN = 2,FUN=sd))
    X.std.p.s[LG.GetInterval(l.obs,is),] <- apply(X=X[LG.GetInterval(l.obs,is),],MARGIN = 2,FUN=standardization)
  }

  Y.std.p.s <- Y
  mean.Y <- NULL
  sd.Y <- NULL
  for (is in 1:s){
    mean.Y <- rbind(mean.Y, apply(X=Y[LG.GetInterval(l.obs,is),],MARGIN = 2,FUN=mean))
    sd.Y <- rbind(sd.Y, apply(X=Y[LG.GetInterval(l.obs,is),],MARGIN = 2,FUN=sd))
    Y.std.p.s[LG.GetInterval(l.obs,is),] <- apply(X=Y[LG.GetInterval(l.obs,is),],MARGIN = 2,FUN=standardization)
  }
  
  
  return(list(X=X.std.p.s,Y=Y.std.p.s,
              mean.X = mean.X, sd.X = sd.X,
              mean.Y = mean.Y, sd.Y = sd.Y))
}

my.perf.spls <- function(keepX,n.train, n.test,l.param,seed,method,method.generation){
  
  # Compute the methods and evaluate their performance (MSEP, TPR, TD) on synthetic data. 
  # Those data are generated following the methodology of the article. 
  
  ### Args
  # keepX : For sPLS, MINT and sPLS for structured data, refers to the number of variables kept. For the others refers to the number of groups kept.
  # n.train :  total number of observations in the training set
  # n.test : total number of observations in the test set
  # l.param : list of parameters for the generation of the synthetic data
  # seed :  seed to reproduce the results
  # method : "spls", "MINT", "sgpls", "spls_structure", "sgpls_structure"
  # method.generation : whether "pb5" for first simulation case and "pb6" for second simualtion case. 
  
  ### Out
  # sum.square : a vector of 4 values giving the MSEP for each vaiable of Y and the overall MSEP
  # TPR : TPR for the selected variables
  # TD : Total Discordance of the the selected variables
  # TP : true positives in terms of  selected variables
  # FN : False negatives in terms of  selected variables
  # FP : False Positives in terms of  selected variables
  # TN : True negatives in terms of  selected variables
  # model : model performed
  
  
  # Setting the parameters
  
  p <- l.param$p
  g <- l.param$g
  s <- l.param$s
  q <- l.param$q
  D <- l.param$D
  
  # Generation of the training data
  
  train.raw <- generate.data.structure(n.obs = n.train,l.param = l.param,seed = seed,method.generation=method.generation,seed.sample.true.C = 2000+seed,D=D) 
  
  list.group.var <- train.raw$list.group.var
  list.group.obs.train <- train.raw$list.group.obs
  
  
  # Standardization of the training data. In spls the standardization is classic, for other methods the standardization is per observation set
  # Setting "l.obs = n.train" makes classic standardization 
  
  if( method == "spls"){
    train <- standardize.by.study(X=train.raw$X,
                                    Y=train.raw$Y,
                                    l.obs = n.train)
  } else {
    train <- standardize.by.study(X=train.raw$X,
                                  Y=train.raw$Y,
                                    l.obs = list.group.obs.train)
  }
  
  colnames(train$X) <- unlist(lapply(1:p, function(x){y<- toString(paste0("X",x));return(y)} ))
  colnames(train$Y) <- unlist(lapply(1:q, function(x){y<- toString(paste0("Y",x));return(y)} ))
  
  # Application of the method
  
  if (method == "sgpls"){
    
      groupX <- subgroupX <- ceiling(1:p / (p/g))
      cv_pls <- list(X = train$X, Y= train$Y, groupX=groupX, subgroupX=subgroupX)
      
      # Grid of parameter among which the tuning is performed.
      sparsities <- cbind(1- (1:9)/10,
                          rep(0,9),
                          (1:9)/10)
      
      model.tune<- tune_sgspls(pls_obj = cv_pls, sparsities = sparsities, group_seq = 4, scale_resp = FALSE,progressBar = FALSE)
      model<- sgspls(X = train$X, Y = train$Y,
                     ncomp = 1,
                     mode = "regression",
                     keepX = keepX, 
                     scale.x = FALSE,
                     groupX = groupX,
                     subgroupX = subgroupX, 
                     indiv_sparsity_x =  model.tune$parameters$indiv_sparsity_x,
                     subgroup_sparsity_x = 0)

    } else if( (method == "spls_structure")||(method == "sgpls_structure")){
    
    # For method "for structure data", equivalence (21) of the article is used
      
    # Translation
    
    X.train.shift <- translation.shift(train$X,s,list.group.obs = list.group.obs.train)
    Y.train.shift <- translation.shift(train$Y,s,list.group.obs = list.group.obs.train)
    
    
    if (method == "spls_structure"){
    groupX<- ceiling(1:(3*p)/3)
    subgroupX<- ceiling(1:(3*p)/3)
        model <- sgspls(X.train.shift , Y.train.shift , ncomp=1, mode = "regression",
                   keepX = 60, groupX = groupX, scale.x = FALSE, subgroupX = subgroupX, indiv_sparsity_x = 0 , subgroup_sparsity_x =0)
    }
    if (method == "sgpls_structure"){
      groupX<- ceiling(1:(3*p)/(3*p/g))
      subgroupX<- ceiling(1:(3*p)/3)
      
      # Grid of parameter among which the tuning is performed.
      sparsities <- cbind(1- (1:9)/10,
                          (1:9)/10,
                          rep(0,9)
                          )
      
      cv_pls <- list(X = X.train.shift, Y = Y.train.shift, groupX=groupX, subgroupX=subgroupX)
      # Tunning of the parameter
      model.tune<- tune_sgspls(pls_obj = cv_pls, sparsities = sparsities, group_seq = 4, scale_resp = FALSE,progressBar = FALSE)
      # Model with tuned parameter
      model<- sgspls(X = X.train.shift, Y = Y.train.shift,
                     ncomp = 1,
                     mode = "regression",
                     keepX = keepX, 
                     scale.x = FALSE,
                     scale.y =  FALSE,
                     groupX = groupX,
                     subgroupX = subgroupX, 
                     indiv_sparsity_x =  0,
                     subgroup_sparsity_x = model.tune$parameters$subgroup_sparsity_x)
      
    }
  } else {
    # Case of spls and MINT
    # The difference between them is that data have been standardized per observation set in MINT
  model <- spls(train$X,
                train$Y,
                ncomp=1,
                mode = "regression"
               ,keepX=keepX
               ,scale = FALSE
               )
  }
  
  # Generation of the test set
  test.raw <- generate.data.structure(n.obs = n.test,l.param = l.param,seed = seed +1000,method.generation=method.generation,seed.sample.true.C = 2000+seed,D=D)
  list.group.obs.test <- test.raw$list.group.obs
  
  # Standardization in the same way than for the training set
  if (method == "spls"){
  test <- standardize.by.study(X=test.raw$X,
                                 Y=test.raw$Y,
                                 l.obs = n.test)
  } else {
  test <- standardize.by.study(X=test.raw$X,
                                  Y=test.raw$Y,
                                  l.obs = list.group.obs.test)
  }
  
  colnames(test$X) <- unlist(lapply(1:p, function(x){y<- toString(paste0("X",x));return(y)} ))
  colnames(test$Y) <- unlist(lapply(1:q, function(x){y<- toString(paste0("Y",x));return(y)} ))
  
  
  if( (method == "spls_structure")||(method == "sgpls_structure")){
    # Compute the equivalent matrices for "for structured data" methods
    X.test.shift <- translation.shift(test$X,s,list.group.obs = list.group.obs.test)
    res <- predict(model, newdata = X.test.shift,scale = FALSE)
  } else {
    res <- predict(model, newdata = test$X,scale = FALSE)
  }
  
  
  # Compare selected variables to true variables in order to compute TPR and TD
  
  true.C <- c(1:15, 21:35, 41:55,61:75)
  
  if( (method == "spls_structure")||(method == "sgpls_structure")){
    selected.C <- unique(ceiling(which(model$weights$X !=0)/s))
  } else if(method == "sgpls"){
    selected.C <- which(model$weights$X !=0)
  } else {
    selected.C <- which(model$loadings$X !=0)
  }


  TP <- length(which(true.C %in% selected.C))
  FN <- length(which(true.C %in% setdiff(1:p,selected.C)))
  FP <-  length(which(setdiff(1:p,true.C)  %in% selected.C))
  TN <-  length(which(setdiff(1:p,true.C)  %in% setdiff(1:p,selected.C)))

  # TPR
  TPR <- TP/(TP+FN)

  # TD

  TD <- (FP + FN)



  # Calculus sum square
  if( (method == "spls_structure")||(method == "sgpls_structure")||(method == "sgpls")){
    Y.pred <- res[,,1]
  } else {
    Y.pred <- res$predict[,,1]
  }
 
  Y.true <- test.raw$Y


  if( (method == "spls_structure")||(method == "sgpls_structure")){
    # reverse the transformation into equivalent matrices if needed
    Y.pred <- un.translation.shift(Y.pred,s,list.group.obs.test)
    Y.pred.unnormed <- Y.pred
    for (iq in 1:q){
      for (is in 1:s){
        Y.pred.unnormed[LG.GetInterval(list.group.obs.test,is),iq] <- (Y.pred[LG.GetInterval(list.group.obs.test,is),iq]*train$sd.Y[is,iq]) + train$mean.Y[is,iq]
      }
    }
  } else if ((method == "MINT")||(method == "sgpls")){
    Y.pred.unnormed <- Y.pred
    for (iq in 1:q){
      for (is in 1:s){
        Y.pred.unnormed[LG.GetInterval(list.group.obs.test,is),iq] <- (Y.pred[LG.GetInterval(list.group.obs.test,is),iq]*train$sd.Y[is,iq]) + train$mean.Y[is,iq]
      }

    }

  } else {
    Y.pred.unnormed <- Y.pred
    for (iq in 1:q){
      Y.pred.unnormed[,iq] <- (Y.pred[,iq]*train$sd.Y[iq]) + train$mean.Y[iq]
    }
  }


  sum.square <- 1:4
  sum.square[1] <- mean((Y.true[,1] - Y.pred.unnormed[,1])^2)
  sum.square[2] <- mean((Y.true[,2] - Y.pred.unnormed[,2])^2)
  sum.square[3] <- mean((Y.true[,3] - Y.pred.unnormed[,3])^2)
  sum.square[4] <- mean((Y.true - Y.pred.unnormed)^2)

    
  return(list(
    sum.square=sum.square,
              TPR =TPR, TD = TD,
              TP = TP, FN = FN,FP = FP,
              TN = TN,
              model = model
              # ,predict=res
              # ,test = test, train = train,
              # Y.pred = Y.pred,
              # Y.pred.unnormed = Y.pred.unnormed,
              # Y.true = Y.true,
              # train.raw = train.raw,
    # model.tune = model.tune,
    # test.raw =test.raw
    ))
  
}
