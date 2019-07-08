#!/usr/bin/env R

# Run 1 iteration of Boruta permutations

borutaboot_1iter_lasso <- function(labeltoken, seed){
  fnstem = "borutadat_lasso_iter" # filename stem
  library(SummarizedExperiment)
  library(Boruta)
  
  message("loading rda objects and functions...")
  load("sesetfilt_degseahack_targetaml.rda")
  
  # Importance Functions
  impLasso <- function(df, classes, seed=2019){
    # df : data frame to parse, rownames = classifier groupings, colnames = feature ids
    require(glmnet)
    require(SummarizedExperiment)
    set.seed(seed) 
    var.classifier <- response <- as.character(classes)
    y <- factor(response); names(y) <- rownames(df) # response var obj
    x = as.matrix(df) # genes of interest
    contrast <- contrasts(y)
    grid <- 10^ seq(10,-2, length=100)
    
    # use cross-validation on the training model.CV only for lambda
    message("performing cross-validation...")
    cv.fit <- cv.glmnet(x, y, family = "binomial",
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
                 getImp = impLasso)
  
  # save data
  message("completed boruta permutations. Saving...")
  save(bdat, file= paste0(fnstem,labeltoken,".rda",collapse=""))
  message("returning...")
}

suppressWarnings(borutaboot_1iter_lasso(labeltoken = commandArgs(T)[1],
                                      seed = commandArgs(T)[2]))