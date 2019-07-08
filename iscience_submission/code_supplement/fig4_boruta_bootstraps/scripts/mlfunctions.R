# ML functions

# Importance Functions
impLasso <- function(df, classes, trainindices, seed=2019){
  # df : data frame to parse, rownames = classifier groupings, colnames = feature ids
  require(glmnet)
  require(SummarizedExperiment)
  set.seed(seed) 
  var.classifier <- response <- as.character(classes)
  y <- factor(response); names(y) <- rownames(df) # response var obj
  x = df # genes of interest
  contrast <- contrasts(y)
  grid <- 10^ seq(10,-2, length=100)
  
  # use cross-validation on the training model.CV only for lambda
  message("performing cross-validation...")
  cv.fit <- cv.glmnet(x[trainindices,], y[trainindices], family = "binomial",
                      type.logistic="modified.Newton", standardize = FALSE,
                      lambda = grid, alpha=1, nfolds = length(trainindices), #LOOCV 
                      type.measure = "class", intercept = FALSE)
  #Select lambda min.
  message("selecting lambda min...")
  lambda.min <- cv.fit$lambda.min
  
  #Fit the full dataset.
  message("fitting whole dataset")
  lassofit <- glmnet(x, y, family = "binomial", standardize = FALSE, 
                     lambda = grid, alpha = 1, intercept = FALSE)
  
  #Extract the coefficients
  lassoimp <- predict(lassofit, type="coefficients", s=lambda.min)
  lassoimp <- lassoimp[2:nrow(lassoimp),1]
  return(lassoimp)
}
impRF <- function(df, classes, ntrees=100, seed=2019){
  require(randomForest)
  set.seed(seed)
  class <- as.numeric(classes)
  rffit <- randomForest(class ~ ., data = as.matrix(df), ntree = ntrees,proximity = TRUE)
  # rfimp <- as.numeric(getRFIvar(rfmodel=rffit))
  rfimp <- importance(rffit)[,1]
  return(rfimp)
}
impXGB <- function(df, classes, seed=2019){
  require(xgboost)
  set.seed(2019)
  message("fitting xgboost model...")
  xgbfit <- xgboost(data = df, label = classes, max_depth = 2,
                    eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic")
  message("getting xgboost importances...")
  xgbimp <- xgb.importance(feature_names = colnames(df), model = xgbfit)
  message("reformatting xgb importances...")
  
  xgbimp.format <- c(rep(0, ncol(df)))
  names(xgbimp.format) <- colnames(df)
  for(f in 1:nrow(xgbimp)){
    xgbimp.format[which(colnames(df)==as.character(xgbimp[f,1]))] <- as.numeric(xgbimp[f,2])
  }
  xgbimp.format <- as.numeric(xgbimp.format)
  names(xgbimp.format) <- colnames(df)
  
  return(xgbimp.format)
}
impSVM <- function(df, classes, seed=2019){
  require(e1071)
  set.seed(2019)
  message("fitting svm model...")
  svmfit <- svm(as.factor(classes)~., 
                data=df, 
                method="C-classification", 
                kernel="linear")
  svmimp <-  t(svmfit$coefs) %*% svmfit$SV
  message("reformatting importances for output...")
  svmimp.format <- svmimp[1,]
  names(svmimp.format) <- colnames(svmimp)
  return(svmimp.format)
}

# Consensus Importance
impCML <- function(df, classes, trainindices, seed, algo.opt=c("lasso","rf","svm","xgb"), 
                   standtable=FALSE, ranksummary="median"){
  # impCML
  # Get consensus importance ranks from disparate algorithms.
  # Arguments
  # * df (matrix) : data table (rows = instances, columns = features)
  # * classes (character): categorizations of instances
  # * algo.opt (list): List of valid algorithms to use for consensus
  # * seed (int): Set the seed for reproducibility
  # * trainindices (numeric, optional): Vector of dataset (rows) corresponding to training sample subsert. Only used for Lasso cross validation step.
  # * ranksummary (string): Either "median" or "mean", the operation used to calculate consensus rank.
  # Returns:
  # * Consensus rank, optionally a standard output table of importances for selected algorithms (if standtable==TRUE) 
  
  
  print("getting importances...")
  implist <- list()
  implabellist <- c()
  if("lasso" %in% algo.opt){
    imp.lasso <- impLasso(df=as.matrix(df),trainindices=trainindices,classes=classes)
    implist[["lasso"]] <- imp.lasso
    implabellist <- c(implabellist, "lasso")
  }
  if("rf" %in% algo.opt){
    imp.rf <- impRF(df=df, classes=classes)
    implist[["rf"]] <- imp.rf
    implabellist <- c(implabellist, "rf")
  }
  if("xgb" %in% algo.opt){
    imp.xgb <- impXGB(df=df, classes=classes)
    implist[["xgb"]] <- imp.xgb
    implabellist <- c(implabellist, "xgb")
  }
  if("svm" %in% algo.opt){
    imp.svm <- impSVM(df=df, classes=classes)
    implist[["svm"]] <- imp.svm
    implabellist <- c(implabellist, "svm")
  }
  # match importance feature label ordering
  if(length(implist)>1){
    for(i in 1:(length(implist)-1)){
      implist[[i]] <- implist[[i]][order(match(names(implist[[i]]),names(implist[[i+1]])))]
    }
  }
  
  # get importance ranks and rankvars
  message("computing rank importances...")
  impranklist <- list()
  for(i in 1:length(implist)){
    # always compute absolute ranks
    impranklist[[names(implist)[i]]] <- rank(abs(implist[[i]]))
  }
  
  medrank <- c()
  meanrank <- c()
  for(r in 1:length(impranklist[[1]])){
    rvalr <- c()
    for(i in 1:length(impranklist)){
      rvalr <- c(rvalr, impranklist[[i]][r])
    }
    medrank <- c(medrank, median(rvalr))
    meanrank <- c(meanrank, mean(rvalr))
  }
  if(ranksummary=="median"){
    lr <- medrank
  } else{
    lr <- meanrank
  }
  
  if(standtable==TRUE){
    # form the stdout table
    message("forming standard table...")
    st <- matrix(nrow=length(impranklist[1]), ncol=0)
    for(r in 1:length(implist)){
      st <- cbind(st, matrix(implist[r], ncol=1))
      colnames(st)[r] <- implabellist[r]
    }
    st <- cbind(st, matrix(medrank, ncol=1)); colnames(st)[ncol(st)] <- "medrank"
    st <- cbind(st, matrix(meanrank, ncol=1)); colnames(st)[ncol(st)] <- "meanrank"
    
    lr <- list("ranksummary"=lr,
               "standtable"=st)
  }
  message("completed all tasks! Returning...")
  return(lr)
}

# CML for Boruta Implementation
impBorutaCML <- function(x=x, y=y, seed, 
                         algo.opt=c("lasso","rf","svm","xgb"), 
                   ranksummary="median"){
  # impCML
  # Get consensus importance ranks from disparate algorithms.
  # Arguments
  # * df (matrix) : Entire dataset matrix (rows = instances, columns = features)
  # * classes (character): categorizations of instances
  # * algo.opt (list): List of valid algorithms to use for consensus
  # * seed (int): Set the seed for reproducibility
  # * ranksummary (string): Either "median" or "mean", the operation used to calculate consensus rank.
  # Returns:
  # * Consensus rank, optionally a standard output table of importances for selected algorithms (if standtable==TRUE) 
  
  # define the test df if provided or if FALSE
  df.test <- x
  message("proceeding with df of size ",nrow(df.test))
  classes <- y
  print("getting importances...")
  implist <- list()
  implabellist <- c()
  if("lasso" %in% algo.opt){
    # define the trainindices if none provided or FALSE
    message("performing lasso...")
    imp.lasso <- impLasso(df=as.matrix(df.test), 
                          trainindices=seq(1,nrow(df.test)), 
                          classes=classes)
    implist[["lasso"]] <- imp.lasso
    implabellist <- c(implabellist, "lasso")
  }
  if("rf" %in% algo.opt){
    imp.rf <- impRF(df=as.matrix(df.test), classes=classes)
    implist[["rf"]] <- imp.rf
    implabellist <- c(implabellist, "rf")
  }
  if("xgb" %in% algo.opt){
    imp.xgb <- impXGB(df=as.matrix(df.test), classes=classes)
    implist[["xgb"]] <- imp.xgb
    implabellist <- c(implabellist, "xgb")
  }
  if("svm" %in% algo.opt){
    imp.svm <- impSVM(df=as.matrix(df.test), classes=classes)
    implist[["svm"]] <- imp.svm
    implabellist <- c(implabellist, "svm")
  }
  # match importance feature label ordering
  if(length(implist)>1){
    for(i in 1:(length(implist)-1)){
      implist[[i]] <- implist[[i]][order(match(names(implist[[i]]),names(implist[[i+1]])))]
    }
  }
  
  # get importance ranks and rankvars
  message("computing rank importances...")
  impranklist <- list()
  for(i in 1:length(implist)){
    # always compute absolute ranks
    impranklist[[names(implist)[i]]] <- rank(abs(implist[[i]]))
  }
  
  medrank <- c()
  meanrank <- c()
  for(r in 1:length(impranklist[[1]])){
    rvalr <- c()
    for(i in 1:length(impranklist)){
      rvalr <- c(rvalr, impranklist[[i]][r])
    }
    medrank <- c(medrank, median(rvalr))
    meanrank <- c(meanrank, mean(rvalr))
  }
  if(ranksummary=="median"){
    lr <- medrank
  } else{
    lr <- meanrank
  }
  message("completed all tasks! Returning...")
  lr
  #return(lr)
}


