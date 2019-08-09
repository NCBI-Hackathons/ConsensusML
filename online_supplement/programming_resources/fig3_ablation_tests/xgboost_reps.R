require(xgboost)
require(SummarizedExperiment)
require(ComplexHeatmap)
require(circlize)
require(reshape)

degfiltset.name <- "sesetfilt_degseahack_targetaml.rda"
load(degfiltset.name)
seset <- degfilt.se
rownames(seset) <- paste0(rownames(seset),"; ",rowData(seset)$hgnc_symbol)

#--------------
# XGBoost reps
#--------------
# Notes: use hyperparameter set from rep3 of hyperparam. optimization
xgb.data <- t(assay(degfilt.se[,degfilt.se$exptset.seahack=="train"]))
xgb.label <- degfilt.se[,degfilt.se$exptset.seahack=="train"]$deg.risk
max.depth = 50
eta = 1
nthread = 2
nrounds = 50
objective = "binary:logistic"
# test
data.test <- t(assay(degfilt.se[,degfilt.se$exptset.seahack=="test"]))
test.label = degfilt.se[,degfilt.se$exptset.seahack=="test"]$deg.risk

#------------------
# Run XGBoost Reps
#------------------
nreps = 100
xgb.fitlist = list()
genes.exclude <- c()
for(i in 1:nreps){
  xgb.data.train <- xgb.data[,!colnames(xgb.data) %in% genes.exclude]
  xgi <- xgboost(data = xgb.data.train, label = xgb.label, 
                 max_depth = max.depth, eta = eta, nthread = nthread, nrounds = nrounds, 
                 objective = objective)
  xgii <- xgb.importance(feature_names = colnames(xgb.data.train), model = xgi)
  xgb.fitlist[[paste0("iter",i)]] <- list(model=xgi, imp=xgii, genes.inc=colnames(xgb.data.train))
  genes.exclude <- c(genes.exclude, xgii$Feature)
  message(i)
}

#-------------------------
# assess rep performances
#-------------------------
# true NLR
postot <- ncol(degfilt.se[,degfilt.se$exptset.seahack=="test" & degfilt.se$deg.risk==1])
# true LR
falsetot <- ncol(degfilt.se[,degfilt.se$exptset.seahack=="test" & degfilt.se$deg.risk==0])
# return list
fiteval <- list()
for(i in 1:length(xgb.fitlist)){
  # assign corresponding gene data subset to match model i
  data.test.filt <- data.test[,colnames(data.test) %in% xgb.fitlist[[i]]$genes.inc]
  # test set predictions
  predi <- predict(xgb.fitlist[[i]]$model, data.test.filt) # testset predictions
  pbi <- ifelse(predi>0.5, 1, 0) # binarized probabilities
  cmi = table(pbi, test.label) # confusion matrix
  # performance metrics
  tp = cmi[2,2]
  fp = cmi[2,1]
  tn = cmi[1,1]
  fn = cmi[1,2]
  mean.err <- mean(as.numeric(predi > 0.5) != test.label) # mean.err
  tpri = (cmi[2,2]/postot)
  tnri = (cmi[1,1]/falsetot)
  fdri = 1-(tp/(tp+fp))
  forvali = 1-(tn/(tn+fn))
  xgb.fitlist[[i]]$performance <- list("confusionMatrix"=cmi,
                                                  "mean_err"=mean.err,
                                                  "tpr"=tpri,
                                                  "tnr"=tnri,
                                                  "fdr"=fdri,
                                                  "for"=forvali)
  message(i)
}

# testset performances df
dfxgp <- matrix(c(unlist(xgb.fitlist[[1]]$performance[2:6])), nrow=1)
colnames(dfxgp) <- names(xgb.fitlist[[1]]$performance[2:6])
for(i in 2:length(xgb.fitlist)){
  dfxgp <- rbind(dfxgp,matrix(xgb.fitlist[[i]]$performance[2:6], nrow=1))
}
dfxgp <- as.data.frame(dfxgp)
colnames(dfxgp)[5] <- "forvar"
rownames(dfxgp) <- c(paste0("rep",seq(1,nrow(dfxgp))))
dfxgp$rep <- seq(1,nrow(dfxgp))

#--------------------------
# graph model performance
#--------------------------

jpeg("modelperf_xgbreps15.jpg", 5,8,units="in",res=400)

par(mfrow=c(2,1),oma=c(2,2,2,1),mar=c(3,3,2,1))

plot(dfxgp$rep, dfxgp$tpr, col="blue", ylim = c(0.6,1), 
     main="XGBoost Fitted Model Performances")
lines(dfxgp$rep, dfxgp$tpr, col="blue")
points(dfxgp$rep, dfxgp$tnr, col="red")
lines(dfxgp$rep, dfxgp$tnr, col="red")
abline(h=0.8,col="black",lty=2,lwd=1)

legend("topright", c("TPR","TNR"), 
       xpd = TRUE, horiz = TRUE,
       bty = "n", pch = c(1,1), lwd=c(1,1),
       col = c("red","blue"), cex = 0.5)

plot(dfxgp$rep, dfxgp$fdr, col="forestgreen", ylim=c(0,0.3))
lines(dfxgp$rep, dfxgp$fdr, col="forestgreen")
points(dfxgp$rep, dfxgp$forvar, col="purple")
lines(dfxgp$rep, dfxgp$forvar, col="purple")
abline(h=0.2,col="black",lty=2,lwd=1)

legend("bottomright", c("FDR","FOR"), 
       xpd = TRUE, horiz = TRUE,
       bty = "n", pch = c(1,1), lwd=c(1,1),
       col = c("forestgreen","purple"), cex = 0.5)

mtext("Rep",side=1,line=3,cex=1)
mtext("Value",side=2,line=1,cex=1, outer=TRUE)
dev.off()

#--------------------
# expr corr heatmaps
#--------------------
# 0. corhm function
compositeCorHM <- function(seset=degfilt.se, corhmname="xgb_hmcorcomp_alldat_rep1234.jpg"){
  cordat <- t(assay(seset))
  cordatsym <- cordat
  colnames(cordatsym) <- paste0(rowData(seset)$ensembl_gene_id,"; ",
                                rowData(seset)$hgnc_symbol)
  cormat.allsamp <- round(cor(cordat, method="spearman"),3)
  cormatsym.allsamp <- round(cor(cordatsym, method="spearman"),3)
  #filt1 <- which(rownames(cormat.allsamp) %in% names(rep1$nonzero.coef))
  #filt2 <- which(colnames(cormat.allsamp) %in% names(rep2$nonzero.coef))
  #filt3 <- which(colnames(cormat.allsamp) %in% names(rep3$nonzero.coef))
  rep1genes <- xgb.fitlist$iter1$imp$Feature
  rep2genes <- xgb.fitlist$iter2$imp$Feature
  rep3genes <- xgb.fitlist$iter3$imp$Feature
  rep4genes <- xgb.fitlist$iter4$imp$Feature
  filt1 <- which(rownames(cormat.allsamp) %in% rep1genes)
  filt2 <- which(rownames(cormat.allsamp) %in% rep2genes)
  filt3 <- which(rownames(cormat.allsamp) %in% rep3genes)
  filt4 <- which(rownames(cormat.allsamp) %in% rep4genes)
  call.21 <- cormatsym.allsamp[filt1,filt2]
  call.31 <- cormatsym.allsamp[filt1,filt3]
  call.41 <- cormatsym.allsamp[filt1,filt4]
  cordat <- t(assay(seset))
  cordatsym <- cordat
  colnames(cordatsym) <- paste0(rowData(seset)$ensembl_gene_id,"; ",
                                rowData(seset)$hgnc_symbol)
  cormat.allsamp <- round(cor(cordat, method="spearman"),3)
  cormatsym.allsamp <- round(cor(cordatsym, method="spearman"),3)
  # colkey
  hm_data = call.21
  breaks=seq(min(hm_data),max(hm_data),0.05)
  hmcol = colorRamp2(breaks,colorRampPalette(c("blue","purple","green","yellow","orange","red"))(n=length(breaks)))
  
  # hm21 -- rep2 genes vs rep1 genes (y-axis)
  hm_data = call.21
  hm21 <- Heatmap(hm_data,col=hmcol,cluster_columns = TRUE,show_heatmap_legend = TRUE,
                  name="Rho", show_row_names = TRUE, show_column_names = TRUE,
                  column_title = "Rep 2 Gene Features", column_dend_reorder = TRUE,
                  row_dend_reorder = TRUE, 
                  heatmap_legend_param = list(color_bar = "continuous"),
                  row_title = "Rep 1 Gene Features",
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8))
  
  # hm31 -- rep3 genes (x-axis) vs rep1 genes (y-axis)
  hm_data = call.31
  hm31 <- Heatmap(hm_data,col=hmcol,cluster_columns = TRUE,show_heatmap_legend = TRUE,
                  name="Rho", show_row_names = TRUE, show_column_names = TRUE,
                  column_title = "Rep 3 Gene Features", column_dend_reorder = TRUE,
                  row_dend_reorder = TRUE, 
                  heatmap_legend_param = list(color_bar = "continuous"),
                  row_title = "Rep 1 Gene Features",
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8))
  
  hm_data = call.41
  hm41 <- Heatmap(hm_data,col=hmcol,cluster_columns = TRUE,show_heatmap_legend = TRUE,
                  name="Rho", show_row_names = TRUE, show_column_names = TRUE,
                  column_title = "Rep 4 Gene Features", column_dend_reorder = TRUE,
                  row_dend_reorder = TRUE, 
                  heatmap_legend_param = list(color_bar = "continuous"),
                  row_title = "Rep 1 Gene Features",
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8))
  
  jpeg(corhmname,15,8,units="in",res=400)
  hm21+hm31+hm41
  dev.off()
}

# 1. all data
compositeCorHM(seset=degfilt.se)