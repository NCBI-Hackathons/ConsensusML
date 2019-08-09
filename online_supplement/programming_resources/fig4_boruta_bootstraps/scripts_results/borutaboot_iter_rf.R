#!/usr/bin/env R

# Run 1 iteration of Boruta permutations

borutaboot_1iter_rf <- function(labeltoken, seed){
  fnstem = "borutadat_rf_iter" # filename stem
  library(SummarizedExperiment)
  library(Boruta)
  
  message("loading rda objects and functions...")
  load("sesetfilt_degseahack_targetaml.rda")

  # Importance Functions
  impRF <- function(df, classes, ntrees=100, seed=2019){
    require(randomForest)
    set.seed(seed)
    class <- as.numeric(classes)
    rffit <- randomForest(class ~ ., data = as.matrix(df), ntree = ntrees,proximity = TRUE)
    # rfimp <- as.numeric(getRFIvar(rfmodel=rffit))
    rfimp <- importance(rffit)[,1]
    return(rfimp)
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
         getImp = impRF)
  
  # save data
  message("completed boruta permutations. Saving...")
  save(bdat, file= paste0(fnstem,labeltoken,".rda",collapse=""))
  message("returning...")
}

suppressWarnings(borutaboot_1iter_rf(labeltoken = commandArgs(T)[1],
                                  seed = commandArgs(T)[2]))