#!/usr/bin/env R

# Run 1 iteration of Boruta permutations

borutaboot_1iter_xgb <- function(labeltoken, seed){
  fnstem = "borutadat_xgb_iter" # filename stem
  library(SummarizedExperiment)
  library(Boruta)
  
  message("loading rda objects and functions...")
  load("sesetfilt_degseahack_targetaml.rda")
  
  # Importance Functions
  impXGB <- function(df, classes, seed=2019){
    require(xgboost)
    set.seed(2019)
    message("fitting xgboost model...")
    xgbfit <- xgboost(data = as.matrix(df), label = classes, max_depth = 2,
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
                 getImp = impXGB)
  
  # save data
  message("completed boruta permutations. Saving...")
  save(bdat, file= paste0(fnstem,labeltoken,".rda",collapse=""))
  message("returning...")
}

suppressWarnings(borutaboot_1iter_xgb(labeltoken = commandArgs(T)[1],
                                        seed = commandArgs(T)[2]))