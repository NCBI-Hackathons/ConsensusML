#!/usr/bin/env R

# Run 1 iteration of Boruta permutations

borutaboot_1iter <- function(labeltoken, seed){
  library(SummarizedExperiment)
  library(Boruta)
  
  message("loading rda objects and functions...")
  load("sesetfilt_degseahack_targetaml.rda")
  load("impfun_lasso.rda")
  load("impfun_rf.rda")
  load("impfun_xgb.rda")
  load("impfun_svm.rda")
  load("borutafun_impcml.rda")
  
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
  # CML for Boruta Implementation
  impBorutaCML <- function(x=x, y=y, seed, ranksummary="median",
                           algo.opt=c("lasso","rf","svm","xgb")){
    # impCML
    # Get consensus importance ranks from disparate algorithms.
    # Arguments
    # * df (matrix) : Entire dataset matrix (rows = instances, columns = features)
    # * classes (character): categorizations of instances
    # * algo.opt (list): List of valid algorithms to use for consensus
    # * seed (int): Set the seed for reproducibility
    # * ranksummary (string): Either "score", "median", or "mean", the operation used to calculate consensus rank.
    # Returns:
    # * Consensus rank, optionally a standard output table of importances for selected algorithms (if standtable==TRUE) 
    
    # define the test df if provided or if FALSE
    df.test <- as.matrix(x)
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
    if(ranksummary %in% c("median","mean")){
      impranklist <- list()
      # for each rank, calculate normalized rank position
      message("computing normalized rank importances...")
      for(i in 1:length(implist)){
        algoname <- names(implist)[i]
        impvals <- implist[[i]]
        # naive rank
        nranki <- rank(abs(impvals))
        # adj1: compress rank lower bounds
        min.ranki <- min(nranki)
        nranki1 <- nranki
        nranki1[nranki1==min.ranki] <- min(nranki1[!nranki1==min.ranki])-1
        # adj2: normalize rank scale (dif.min/range)
        nranki2 <- nranki1
        nranki2 <- (nranki2-min(nranki2))/(max(nranki2)-min(nranki2))
        # always compute absolute ranks
        impranklist[[algoname]] <- nranki2
        # plot(impvals, nranki, main=algoname)
        # plot(impvals, nranki1, main=algoname)
        # plot(impvals, nranki2, main=algoname)
      }
      medrank <- c()
      meanrank <- c()
      # iterate on features
      for(r in 1:ncol(x)){
        rvalr <- c() # get feature-level ranks
        for(i in 1:length(impranklist)){
          rvalr <- c(rvalr, impranklist[[i]][r]) # append ranks
        }
        # append summarized feature-level ranks, to return
        medrank <- c(medrank, median(rvalr))
        meanrank <- c(meanrank, mean(rvalr))
      }
      # final importance value is directly proportional to summarized rank
      if(ranksummary=="median"){
        lr <- medrank
      } else{
        lr <- meanrank
      }
      message("completed all tasks! Returning...")
      lr
      return(lr)
    }
    
    # get score importance metric
    if(ranksummary=="score"){
      # notes: by-algo score criteria
      # lasso, coeff!=0
      # rf, imp>=90th quantile
      # svm, abs-weight >= 90th quantile
      # xgb, coeff != 0
      featnames <- colnames(x)
      impscore <- rep(0,length(featnames))
      names(impscore) <- colnames(x)
      for(i in 1:length(implist)){
        ni <- names(implist)[i]
        if(ni=="lasso"){
          which.features <- names(implist[[i]][implist[[i]]!=0])
          impscore[which.features] <- impscore[which.features]+1
        }
        if(ni=="xgb"){
          which.features <- names(implist[[i]][implist[[i]]!=0])
          impscore[which.features] <- impscore[which.features]+1
        }
        if(ni=="rf"){
          qrf <- as.numeric(quantile(as.numeric(implist[[i]]), seq(0,1,0.1))[10])
          which.features <- names(implist[[i]][implist[[i]]>=qrf])
          impscore[which.features] <- impscore[which.features]+1
        }
        if(ni=="svm"){
          qsvm <- as.numeric(quantile(abs(as.numeric(implist[[i]])), seq(0,1,0.1))[10])
          which.features <- names(implist[[i]][abs(implist[[i]])>=qsvm])
          impscore[which.features] <- impscore[which.features]+1
        }
      }
      rs <- as.numeric(impscore)
      rs
      return(rs)
    }
    
    return()
  }
  
  message("extracting data and classes...")
  classes <- as.character(degfilt.se$deg.risk)
  data <- t(assay(degfilt.se))
  it <- list() # iteration list with by-gene summaries
  si <-  seq(1, nrow(data), 1) # sample (row) indices
  si.rg0 <- si[which(degfilt.se$deg.risk==0)]
  si.rg1 <- si[which(degfilt.se$deg.risk==1)]
  
  # get random sample subset
  message("retrieving randomized data subset objects...")
  set.seed(seed)
  train.rg0 <- sample(si.rg0, 40)
  train.rg1 <- sample(si.rg1, 51)
  trainindices <- c(train.rg0, train.rg1)
  dfi <- data[trainindices,]
  classesi <- classes[trainindices]
  
  # run boruta permutations
  message("running boruta permutations...")
  bdat <- Boruta(x = dfi, 
         y = classesi, 
         getImp = impBorutaCML)
  
  # save data
  message("completed boruta permutations. Saving...")
  save(bdat, file= paste0("borutadat_",labeltoken,".rda",collapse=""))
  message("returning...")
}

suppressWarnings(borutaboot_1iter(labeltoken = commandArgs(T)[1],
                                  seed = commandArgs(T)[2]))