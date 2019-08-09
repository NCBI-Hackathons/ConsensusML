# make composite importance plot for main figure 2
library(SummarizedExperiment)

load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/lasso_resultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/svm4reps_resultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/rf_noboost_2k5k10ktrees_allresultslist.rda")
load("~/scratch/consensusML/manuscript_draft/programming_resources/hyperparam_optimization/data/xgb_resultslist.rda")

gene.lab <- data.frame(rowData(degfilt.se))
eidvar = rownames(gene.lab)
gene.lab <- data.frame(eid=rownames(gene.lab),
                       symbol=gene.lab$hgnc_symbol,
                       clab=paste0(rownames(gene.lab),":",gene.lab$hgnc_symbol),
                       stringsAsFactors = F)
rownames(gene.lab) <- eidvar

pfheight = 4; pfwidth = 10

# lasso
lr <- lasso.resultslist
lr[[1]]$nonzero.coef
pdf("lasso-betaimp_3reps.pdf",pfwidth,pfheight)
par(mfrow=c(1,3), oma = c(10,2,1,1))
for(i in 1:3){
  bci = lr[[i]]$nonzero.coef; bci <- bci[order(bci)]
  labi = gene.lab[names(bci),]$clab
  barplot(bci, names=labi, col=ifelse(bci<0, "blue", "red"),
          las=2, cex.names = 0.9, main=paste0("Lasso Rep ",i),
          ylab="Importance (Coeff.)")
}
dev.off()

# random forest
rfl = rf.returnlist
cv = c("orange","green","purple")
pdf("rfimp_3reps.pdf", pfwidth, pfheight)
par(mfrow=c(1,3), oma = c(10,2,1,1))
for(i in 1:3){
  rfi = rfl[[i]]$fitmodel$importance
  qfi = quantile(rfi[,1], seq(0,1,0.01))[100] # 99%
  rfif = rfi[rfi[,1]>=qfi,]; rfif <- rfif[order(rfif)]
  labi = gene.lab[names(rfif),]$clab
  barplot(rfif, 
          names=labi, 
          col=cv[i],
          las=2, cex.names = 0.9, 
          main=paste0("RF Rep ",i),
          ylab="Importance (Mean Dev. Gini)")
}
dev.off()

# svm
svl = svm.resultslist
#cv = c("orange","green","purple","pink")
pdf("svmimp_3reps.pdf", pfwidth, pfheight)
par(mfrow=c(1,4), oma = c(10,2,1,1))

for(i in 1:4){
  swi = t(svl[[i]]$weightsvect)
  qfi = quantile(abs(swi[,1]), seq(0,1,0.01))[100] # 99% abs
  swif = swi[,1]; names(swif) <- rownames(swi)
  swif = swif[swif>=qfi]; 
  swif <- swif[order(swif)]
  labi = gene.lab[names(swif),]$clab
  barplot(swif, 
          names=labi, 
          col=ifelse(swif>0,"red","blue"),
          las=2, cex.names = 0.9, 
          main=paste0("SVM Rep ",i),
          ylab="Importance (Weight)")
}
dev.off()

# xgboost
xgl = xg.resultslist[1:3]
cv = c("orange","green","purple","pink")
pdf("xgbimp_3reps.pdf", pfwidth, pfheight)
par(mfrow=c(1,3), oma = c(10,2,1,1))
for(i in 1:3){
  xgi = as.data.frame(xgl[[i]]$importance)
  xgi2 <- as.numeric(xgi[,"Gain"])
  names(xgi2) <- xgi$Feature
  xgi2 <- xgi2[order(xgi2)]
  barplot(xgi2, 
          names=names(xgi2), 
          col=cv[i],
          las=2, cex.names = 0.9, 
          main=paste0("XGB Rep ",i),
          ylab="Importance (Gain)")
}
dev.off()
