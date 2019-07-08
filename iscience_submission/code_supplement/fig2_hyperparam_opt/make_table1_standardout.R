# make the standard output table, table 1

library(SummarizedExperiment)

load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/sesetfilt_degseahack_targetaml.rda")

load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/lasso_resultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/svm4reps_resultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/rf_noboost_2k5k10ktrees_allresultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/xgb_resultslist.rda")

# make standardized output table
st <- as.data.frame(rowData(degfilt.se))
st$ensembl_id <- rownames(st)
st <- st[,c("ensembl_id","hgnc_symbol")]

imp.s <- rep(0,nrow(st))
names(imp.s) <- st$ensembl_id

# append lasso results
for(i in 1:3){
  li = lasso.resultslist[[i]]
  lic <- li$nonzero.coef
  li.imp <- imp.s
  for(k in 1:length(lic)){
    li.imp[names(li.imp)==names(lic)[k]] <- lic[k]
  }
  st$newcol <- li.imp
  colnames(st)[ncol(st)] <- paste0("lasso_rep",i)
}

# append rf results
for(i in 1:3){
  rfii <- rf.returnlist[[i]]$fitmodel$importance
  rf.imp <- imp.s
  for(k in 1:nrow(rfii)){
    rf.imp[rownames(rfii)[k]] <- rfii[k,1]
  }
  st$newcol <- rf.imp
  colnames(st)[ncol(st)] <- paste0("rf_rep",i)
}

# append xgb results
for(i in 1:3){
  xgii <- as.data.frame(xg.resultslist[[i]]$importance)
  xg.imp <- imp.s
  for(k in 1:nrow(xgii)){
    xg.imp[gsub(";.*","",as.character(xgii[k,1]))] <- xgii[k,2]
  }
  st$newcol <- xg.imp
  colnames(st)[ncol(st)] <- paste0("xgb_rep",i)
}

# append svm results
for(i in 1:4){
  svmii <- t(svm.resultslist[[i]]$weightsvect)
  svm.imp <- imp.s
  for(k in 1:nrow(svmii)){
    svm.imp[rownames(svmii)[k]] <- svmii[k,1]
  }
  st$newcol <- svm.imp
  colnames(st)[ncol(st)] <- paste0("svm_rep",i)
}

write.csv(st, "table2_standardoutput.csv")
save(st, file="stout_table.rda")


