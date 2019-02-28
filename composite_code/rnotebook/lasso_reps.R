# lasso reps and model assessment
# Note: reps of lasso were performed on DEGs for TARGET pediatric AML samples
#   Models were optimized using the training subset, then assessed using the validation subset
#   Lasso iterations describe runs after performing lasso on subsets of initial DEGs (>1st rep)

# define globals
sys.sep = "/"
data.dir = "data"
seobj.dir = "seobjects"
figs.dir = "figures"

# standard outputs table
stdtable.name <- "standardtable_mloutputs_summary.rda"

# Load the primary summarized experiment object for experiment
degfiltset.name <- "sesetfilt_degseahack_targetaml.rda"
load(paste0(data.dir, sys.sep, degfiltset.name))

# lasso function
runLasso <- function(seset, seed=2019){
  # runLasso
  # Fit a model using penalized regression with lasso
  # Arguments:
  # * sese: Valid summarized experiment object
  # * seed: (int) set seed for randomization
  # Returns:
  # * resultslist (list) : Results of lasso fit
  require(glmnet)
  require(SummarizedExperiment)
  set.seed(seed) 
  
  gene.names = as.character(rownames(rowData(seset)))
  var.classifier = seset$deg.risk
  df = t(assay(seset))
  train.names = colnames(assay(seset[,seset$exptset.seahack=="train"]))
  test.names = colnames(assay(seset[,seset$exptset.seahack=="test"]))
  response <- var.classifier
  predictors <- gene.names
  y <- factor(response); names(y) <- colnames(assay(seset)) # response var obj
  x = df[,colnames(df) %in% predictors] # genes of interest
  contrast <- contrasts(y)
  grid <- 10^ seq(10,-2, length=100)
  standardize = FALSE
  fit <- glmnet(x[train.names,], y[train.names], family = "binomial", alpha=1, 
                standardize = standardize, lambda = grid, intercept = FALSE)
  
  # use cross-validation on the training model.CV only for lambda
  cv.fit <- cv.glmnet(x[train.names,], y[train.names], family = "binomial",
                      type.logistic="modified.Newton", standardize = standardize,
                      lambda = grid, alpha=1, nfolds = length(train.names), #LOOCV 
                      type.measure = "class", intercept = FALSE)
  #Select lambda min.
  lambda.min <- cv.fit$lambda.min
  #predict the classes
  pred.class <- predict(fit, newx = x[test.names,], type="class", s=lambda.min)
  #find the test error
  tab <- table(pred.class,y[test.names])
  testError <- mean(pred.class != y[test.names]) #how many predicted classes were incorrect
  #Fit the full dataset.
  final <- glmnet(x, y,family = "binomial", standardize = standardize, 
                  lambda = grid, alpha = 1, intercept = FALSE)
  #Extract the coefficients
  coef <- predict(final, type="coefficients", s=lambda.min)
  idx <- which(!as.numeric(coef)==0)
  nonZero <- coef[idx,]
  #Results 
  resultslist <- list(train.names, test.names, contrast, fit, 
                       cv.fit, tab, testError, final, nonZero, seed)
  names(resultslist) <- c("training.set", "testing.set","contrast", "train.fit",
                           "cv.fit", "confusionMatrix","test.error", "final.model", 
                           "nonzero.coef", "seed")
  return(resultslist)
}

# lasso reps
genes.exclude <- c() # running list of genes (selected features) to exclude

rep1 <- runLasso(seset=degfilt.se)
genes.exclude <- c(genes.exclude, names(rep1$nonzero.coef))

rep2 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep2$nonzero.coef))

rep3 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep3$nonzero.coef))

rep4 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep4$nonzero.coef))

rep5 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep5$nonzero.coef))

rep6 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep6$nonzero.coef))

rep7 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude])
genes.exclude <- c(genes.exclude, names(rep7$nonzero.coef))

rep8 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep8$nonzero.coef))

rep9 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep9$nonzero.coef))

rep10 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep10$nonzero.coef))

rep11 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep11$nonzero.coef))

rep12 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep12$nonzero.coef))

rep13 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep13$nonzero.coef))

rep14 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep14$nonzero.coef))

rep15 <- runLasso(seset=degfilt.se[!rownames(degfilt.se) %in% genes.exclude,])
genes.exclude <- c(genes.exclude, names(rep15$nonzero.coef))

replist <- list(rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8,
                rep9, rep10, rep11, rep12, rep13, rep14, rep15)

# graph model performance

dfp <- as.data.frame(matrix(nrow=15,ncol=5))
colnames(dfp) <- c("rep","tpr","tnr","fdr","for")
dfp$rep <- seq(1,15,1)

postot <- ncol(degfilt.se[,degfilt.se$exptset.seahack=="test" & degfilt.se$deg.risk==1])
falsetot <- ncol(degfilt.se[,degfilt.se$exptset.seahack=="test" & degfilt.se$deg.risk==0])

# True positive rate, TPR = TP/P
tprvar <- c()
for(i in 1:length(replist)){
  tprvar <- c(tprvar, (replist[[i]]$confusionMatrix[2,2]/postot))
}
dfp$tpr <- tprvar

# True negative rate, TNR = TF/F
tnrvar <- c()
for(i in 1:length(replist)){
  tnrvar <- c(tnrvar,(replist[[i]]$confusionMatrix[1,1]/falsetot))
}
dfp$tnr <- tnrvar

# False discovery rate, FDR = 1 - (TP/[TP+FP])
fdrvar <- c()
for(i in 1:length(replist)){
  tp = replist[[i]]$confusionMatrix[2,2]
  fp = replist[[i]]$confusionMatrix[2,1]
  fdrvar <- c(fdrvar, 1-(tp/(tp+fp)))
}
dfp$fdr <- fdrvar

# False omission rate, FOR = 1 - (TN/[TN+FN])
forvar <- c()
for(i in 1:length(replist)){
  tn = replist[[i]]$confusionMatrix[1,1]
  fn = replist[[i]]$confusionMatrix[1,2]
  forvar <- c(forvar, 1-(tn/(tn+fn)))
}
dfp$forvar <- forvar

#=========================
# graph model performance
#=========================

jpeg(paste0(figs.dir, sys.sep, "modelperf_lassoreps15.jpg"), 5,8,units="in",res=400)
par(mfrow=c(2,1),oma=c(2,2,2,1),mar=c(3,3,2,1))

plot(dfp$rep, dfp$tpr, col="blue", ylim = c(0.6,1), 
     main="Lasso Fitted Model Performances")
lines(dfp$rep, dfp$tpr, col="blue")
points(dfp$rep, dfp$tnr, col="red")
lines(dfp$rep, dfp$tnr, col="red")
abline(h=0.8,col="black",lty=2,lwd=1)

legend("bottomleft", c("TPR","TNR"), 
       xpd = TRUE, horiz = TRUE,
       bty = "n", pch = c(1,1), lwd=c(1,1),
       col = c("red","blue"), cex = 1)

plot(dfp$rep, dfp$fdr, col="forestgreen", ylim=c(0,0.3))
lines(dfp$rep,dfp$fdr, col="forestgreen")
points(dfp$rep,dfp$forvar, col="purple")
lines(dfp$rep,dfp$forvar, col="purple")
abline(h=0.2,col="black",lty=2,lwd=1)

legend("topleft", c("FDR","FOR"), 
       xpd = TRUE, horiz = TRUE,
       bty = "n", pch = c(1,1), lwd=c(1,1),
       col = c("forestgreen","purple"), cex = 1)

mtext("Rep",side=1,line=3,cex=1)
mtext("Value",side=2,line=1,cex=1, outer=TRUE)
dev.off()

#========================
# first 3 reps, heatmaps
#========================
require(ComplexHeatmap)
require(circlize)
require(reshape)

#==================
# all samples corr
#==================
# 0. corhm function
compositeCorHM <- function(seset, corhmname="lasso_hmcorcomp_alldat_rep123.jpg"){
  cordat <- t(assay(seset))
  cordatsym <- cordat
  colnames(cordatsym) <- paste0(rowData(seset)$ensembl_gene_id,"; ",
                                rowData(seset)$hgnc_symbol)
  
  cormat.allsamp <- round(cor(cordat, method="spearman"),3)
  cormatsym.allsamp <- round(cor(cordatsym, method="spearman"),3)
  
  filt1 <- which(rownames(cormat.allsamp) %in% names(rep1$nonzero.coef))
  filt2 <- which(colnames(cormat.allsamp) %in% names(rep2$nonzero.coef))
  filt3 <- which(colnames(cormat.allsamp) %in% names(rep3$nonzero.coef))
  
  call.21 <- cormatsym.allsamp[filt1,filt2]
  call.31 <- cormatsym.allsamp[filt1,filt3]
  
  cordat <- t(assay(seset))
  cordatsym <- cordat
  colnames(cordatsym) <- paste0(rowData(seset)$ensembl_gene_id,"; ",
                                rowData(seset)$hgnc_symbol)
  
  cormat.allsamp <- round(cor(cordat, method="spearman"),3)
  cormatsym.allsamp <- round(cor(cordatsym, method="spearman"),3)
  
  filt1 <- which(rownames(cormat.allsamp) %in% names(rep1$nonzero.coef))
  filt2 <- which(colnames(cormat.allsamp) %in% names(rep2$nonzero.coef))
  filt3 <- which(colnames(cormat.allsamp) %in% names(rep3$nonzero.coef))
  
  call.21 <- cormatsym.allsamp[filt1,filt2]
  call.31 <- cormatsym.allsamp[filt1,filt3]
  
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
  
  jpeg(paste0(figs.dir, sys.sep, corhmname),12,8,units="in",res=400)
  hm21+hm31
  dev.off()
}

# 1. all data
compositeCorHM(seset=degfilt.se)

# 2. data subsets
# test subset
compositeCorHM(seset=degfilt.se[,degfilt.se$exptset.seahack=="test"],
               corhmname="lasso_hmcorcomp_alltest_rep123.jpg")
# train subset
compositeCorHM(seset=degfilt.se[,degfilt.se$exptset.seahack=="train"],
               corhmname="lasso_hmcorcomp_alltrain_rep123.jpg")

# 3. all data, risk group subsets
compositeCorHM(seset=degfilt.se[,degfilt.se$deg.risk==0],
               corhmname="lasso_hmcorcomp_alldat0_rep123.jpg")
compositeCorHM(seset=degfilt.se[,degfilt.se$deg.risk==1],
               corhmname="lasso_hmcorcomp_alldat1_rep123.jpg")
