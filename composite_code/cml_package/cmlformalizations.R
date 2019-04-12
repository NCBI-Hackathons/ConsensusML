# CML formalizations for TARGET AML dataset
# Test formalization of consensus approach, and compare to standard approach with Boruta.

library(Boruta)
library(SummarizedExperiment)

load("sesetfilt_degseahack_targetaml.rda")
amldf <- as.data.frame(t(assay(degfilt.se)))
dim(amldf)
# [1]  137 1984
# ownames(amldf) <- degfilt.se$deg.risk
amldf$brg <- as.factor(degfilt.se$deg.risk)

#----------------------------
# Testing Boruta Hyperparams
#----------------------------
# rf importance
ba.rf <- Boruta(brg~., data=amldf)

# other questions:
# does boruta improve on xgboost?
# does boruta improve feature overlap, rf vs xgboost?

#---------------------------------------------------
# Compare CML to RF based importance in Boruta runs
#---------------------------------------------------
df <- t(assay(degfilt.se))
classes <- as.character(degfilt.se$deg.risk)
trainindices = which(degfilt.se$exptset.seahack=="train")

ba.cml <- Boruta(x=as.matrix(df), y = classes, getImp = impBorutaCML)
ba.rf <- Boruta(x=df, y = as.factor(classes), getImp = getImpRfZ)

# table(ba.cml$finalDecision)
# Tentative Confirmed  Rejected 
#  0        12      1972

table(ba.rf$finalDecision)
# Tentative Confirmed  Rejected 
#   89        55      1840

final.cml <- ba.cml$finalDecision
final.rf <- ba.rf$finalDecision

# Consenus of two Baruta trials
length(intersect(names(final.cml[final.cml=="Confirmed"]),
                 names(final.rf[final.rf=="Confirmed"])))
# 12

st <- read.csv("standouttable.csv")
length(intersect(st[!st$lasso_coef_rep1==0,]$X, names(final.cml[final.cml=="Confirmed"])))
# [1] 4
gx <- intersect(st[!st$lasso_coef_rep1==0,]$X, names(final.cml[final.cml=="Confirmed"]))
# "ENSG00000260182.1"  "ENSG00000132975.7"  "ENSG00000078399.14" "ENSG00000267453.5" 
st[st$X %in% c(gx),]

bcml.genes <- names(final.cml[final.cml=="Confirmed"])
st[st$X %in% bcml.genes,]












