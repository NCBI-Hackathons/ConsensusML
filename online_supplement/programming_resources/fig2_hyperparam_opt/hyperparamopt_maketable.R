# Make table summarizing hyperparameter optimization results
library(SummarizedExperiment)

# preamble and load objects
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/sesetfilt_degseahack_targetaml.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/svm4reps_resultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/lasso_resultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/rf_noboost_2k5k10ktrees_allresultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/xgb_resultslist.rda")

svml <- svm.resultslist
lassol <- lasso.resultslist
rfl <- rf.returnlist
xgbl <- xg.resultslist

train.classes <- degfilt.se[,degfilt.se$exptset.seahack=="train"]$deg.risk
test.classes <- degfilt.se[,degfilt.se$exptset.seahack=="test"]$deg.risk

#------------------------------
# retrieve performance metrics
#------------------------------
svm1.cm <- table(svml$svm1$predictions_test, test.classes)
svm2.cm <- table(svml$svm2$predictions_test, test.classes)
svm3.cm <- table(svml$svm3$predictions_test, test.classes)
svm4.cm <- table(svml$svm4$predictions_test, test.classes)

lasso1.cm <- lassol[[1]]$confusionMatrix
lasso2.cm <- lassol[[2]]$confusionMatrix
lasso3.cm <- lassol[[3]]$confusionMatrix
lasso4.cm <- lassol[[4]]$confusionMatrix

rfl1.cm <- t(rfl$rf2k.results$conf.matrix)
rfl2.cm <- t(rfl$rf5k.results$conf.matrix)
rfl3.cm <- t(rfl$rf10k.results$conf.matrix)

xgb1.cm <- xgbl$rep1$performance_testset$confusionMatrix
xgb2.cm <- xgbl$rep2$performance_testset$confusionMatrix
xgb3.cm <- xgbl$rep3$performance_testset$confusionMatrix
xgb4.cm <- xgbl$rep4$performance_testset$confusionMatrix

# note all CM have rows as pred, col as classes
cml <- list(svm=list(i1=svm1.cm,
                     i2=svm2.cm,
                     i3=svm3.cm,
                     i4=svm4.cm),
            lasso=list(i1=lasso1.cm,
                       i2=lasso2.cm,
                       i3=lasso3.cm,
                       i4=lasso4.cm),
            rf=list(i1=rfl1.cm,
                    i2=rfl2.cm,
                    i3=rfl3.cm),
            xgb=list(i1=xgb1.cm,
                     i2=xgb2.cm,
                     i3=xgb3.cm,
                     i4=xgb4.cm))

#------------------------------
# aggregate performance metrics
#------------------------------
cmt <- matrix(nrow=0,ncol=10)
for(ll in 1:length(cml)){
  ll.i <- cml[[ll]]
  namel <- names(cml)[ll]
  cmidat <- c()
  for(i in 1:length(ll.i)){
    cmidat <- c(namel, as.character(i))
    cmi = ll.i[[i]]
    # tn, fn, tp, fp
    cmidat <- c(cmidat, 
                cmi[1,1], cmi[2,1],
                cmi[2,2], cmi[1,2])
    # tpr, tnr, fdr, for 
    cmidat <- c(cmidat, 
                cmi[2,2]/length(test.classes[test.classes==1]), # tpr
                cmi[1,1]/length(test.classes[test.classes==0]), # tnr
                cmi[2,1]/(cmi[2,1]+length(test.classes[test.classes==1])), # fdr
                cmi[1,2]/(cmi[1,2]+length(test.classes[test.classes==0])) # for
                )
    cmt <- rbind(cmt, cmidat)
    message(ll," ",i)
  }
  message(ll)
}
colnames(cmt) <- c("algorithm","rep","TN","FN","TP","FP",
                   "TPR","TNR","FDR","FOR")

