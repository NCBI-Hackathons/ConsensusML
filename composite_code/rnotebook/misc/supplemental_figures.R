# Code for supplemental figures to ConsensusML+Target Pediatric AML RNAseq biomarkers project

# libraries
library(plyr)
library(mlr)
library(magrittr)
library(ggplot2)
library(EnsDb.Hsapiens.v75)
library(glmnet)
library(ROSE)
library(knitr)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(limma)
library(edgeR)
library(MLSeq)
library(DESeq2)
library(xlsx)
library(VennDiagram)
library(SummarizedExperiment)
library(GenomicRanges)
library(circlize)
library(reshape2)

# define globals
sys.sep = "/"
data.dir = "data"
seobj.dir = "seobjects"
figs.dir = "figures"
stdtable.name <- "standardtable_mloutputs_summary.rda"
# degseset.name <- "seset_degseahack_targetaml.rda"
degfiltset.name <- "sesetfilt_degseahack_targetaml.rda"

load(paste0(data.dir, sys.sep, degfiltset.name))
load(paste0(data.dir, sys.sep, stdtable.name))

#=================
# lasso corplots
# Make cor plot of top features selected by 
# initial lasso, and lasso subset 
# (excluding initial selected features)

# corr data, full dataset
cormat <- round(cor(t(assay(degfilt.se)), method="spearman"),3)
melted_cormat <- melt(cormat)
melted_cormat <- melted_cormat[!melted_cormat$Var1==melted_cormat$Var2,] # pre-filter

# corr table, filtered genes
goi <- rownames(standtable[!standtable$lasso_coeff==0 | !standtable$lasso_rmCoeff==0 | !standtable$lasso_rep3_cofilt2==0,])
cormat_filt <- round(cor(t(assay(degfilt.se[goi,])), method="spearman"),3)
melted_cormat_filt <- melt(cormat_filt)
melted_cormat_filt <- melted_cormat_filt[!melted_cormat_filt$Var1==melted_cormat_filt$Var2,]

#====================
# Density Histograms
#====================
corhist.name <- "corplot_hist_deg-vs-lasso3reps.jpg"
jpeg(paste0(figs.dir, sys.sep, corhist.name), 6,4,units="in",res=400)
plot(density(melted_cormat$value), col=rgb(0.2,0.2,0.8,alpha=0.2), lwd=2, lty=2,
     xlim=c(-1.2,1.2), ylab="Relative Density", xlab="Rho", main="Gene Expr. Correlation Dist.")
lines(density(melted_cormat_filt$value), col=rgb(0.1,0.9,0.2),lwd=3, lty=1)
legend("topright",bty="n",legend=c("All DEGs", "Lasso Genes\n(3 reps, with and\nwithout coeff. filt)"), col=c(rgb(0.2,0.2,0.8,alpha=0.2), rgb(0.1,0.9,0.2)), lwd=c(2,3),lty=c(2,1),cex=0.5)
dev.off()

#=======================
# Correlation ggHeatmap
#=======================
# color scale
colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(42)

mcf <- melted_cormat

# rep1 vs rep2 genes
goi1 <- rownames(standtable[!standtable$lasso_coeff==0,])
goi2 <- rownames(standtable[!standtable$lasso_rmCoeff==0,])
mcf_12 <- mcf[mcf$Var2 %in% goi1 & mcf$Var1 %in% goi2,]
colnames(mcf_12) <- c("Rep2 Feature","Rep1 Feature","Rho")
corhm12 = ggplot(data = mcf_12, aes(x=mcf_12$`Rep2 Feature`, 
                                  y=mcf_12$`Rep1 Feature`, 
                                  fill=mcf_12$Rho)) + 
  geom_tile() +
  scale_fill_gradientn(colours = colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Spearman Corr. Matrix") +
  xlab("Rep2 Features") +
  ylab("Rep1 Features") +
  labs(fill="Rho")

# rep1 vs rep3 genes
goi1 <- rownames(standtable[!standtable$lasso_coeff==0,])
goi3 <- rownames(standtable[!standtable$lasso_rep3_cofilt2==0,])
# mcf_13 <- mcf[mcf$Var1 %in% goi1 & mcf$Var2 %in% goi3,]
mcf_13 <- mcf[mcf$Var1 %in% goi1 & mcf$Var2 %in% goi3,]
colnames(mcf_13) <- c("Rep1 Feature","Rep3 Feature","Rho")
corhm13 = ggplot(data = mcf_13, aes(x=mcf_13$`Rep3 Feature`, 
                                    y=mcf_13$`Rep1 Feature`, 
                                    fill=mcf_13$Rho)) + 
  geom_tile() +
  scale_fill_gradientn(colours = colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Rep 3 Features") +
  xlab("Rep3 Features") +
  ylab("Rep1 Features") +
  labs(fill="Rho")

# save plots
corhm.name.12 = "corplot_lasso-rep1-vs-rep2.jpg"
jpeg(paste0(figs.dir, sys.sep, corhm.name.12), 7, 7, units="in", res=400)
corhm12 + scale_colour_gradientn(colours=rainbow(4))
dev.off()

corhm.name.13 = "corplot_lasso-rep1-vs-rep3.jpg"
jpeg(paste0(figs.dir, sys.sep, corhm.name.13), 7, 7, units="in", res=400)
corhm13 + scale_colour_gradientn(colours=rainbow(4))
dev.off()

#================================
# complex heatmap aggregate plot
#================================

# corhm12 + corhm13


