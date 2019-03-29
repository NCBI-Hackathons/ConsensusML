library(SummarizedExperiment)
library(gridExtra)
library(FactoMineR)
library(factoextra)
library(ggbiplot)
library(enrichR)
library(UpSetR)
library(ggcorrplot)
library(Rtsne)
#library(Rtsne)

# Notes on the following analysis:
# 1. "Important" features thresholds in SVM and RF:
# 1A. SVM: feature abs weight >0.008, results in 
# 1B. RF: feature importnce >0.1

load("sesetfilt_degseahack_targetaml.rda")
data = assay(degfilt.se)
st <- read.csv("standouttable.csv")

#------------------------------------
# Principle Component Analysis (PCA)
#------------------------------------
data_pca <- prcomp(t(data), scale.=T)
res.pca <- PCA(t(data), graph = FALSE)

# get gene contributions to components 1 and 1:3
# note: the following is informed by function 'fviz_contrib'
X = data_pca
choice = "var"
axes = 1 # first pc
dd1 <- facto_summarize(X, element = choice, result = "contrib", 
                       axes = axes)
axes = c(1,2,3) # first 3 pcs
dd123 <- facto_summarize(X, element = choice, result = "contrib", 
                         axes = axes)
identical(dd1$name,dd123$name)
colnames(dd1) <- c("name","cont_pca1")
pccont <- cbind(dd1,data.frame(cont_pca123=dd123$contrib))

# append whether feature is above/below uniformity threshold
contrib <- pccont$cont_pca1
names(contrib) <- rownames(dd)
theo_contrib1 <- 100/length(contrib)
axes = c(1,2,3) # first 3 pcs
eig <- get_eigenvalue(X)[axes, 1]
theo_contrib123 <- sum(theo_contrib1 * eig)/sum(eig)
pccont$above_theoimp_pc1 <- ifelse(pccont$cont_pca1>theo_contrib1,"Y","N")
pccont$above_theoimp_pc123 <- ifelse(pccont$cont_pca123>theo_contrib123,"Y","N")
table(pccont$above_theoimp_pc1, pccont$above_theoimp_pc123)
#    N   Y
# N 807 303
# Y 239 635

save(pccont, file="pca_degcont.rda")
write.csv(pccont, file="pca_degcont.csv")

#-----------------------
# Visualize PCA Results
#-----------------------
# pca scatterplots, dims 1 and 2
p1 <- ggbiplot(data_pca, ellipse=TRUE, groups=degfilt.se$Risk.group, var.axes=FALSE) + ggtitle("Risk Group")
p2 <- ggbiplot(data_pca, ellipse=TRUE, groups=degfilt.se$deg.risk, var.axes=FALSE) + ggtitle("Risk Group (binary)")
p3 <- ggbiplot(data_pca, ellipse=TRUE, groups=degfilt.se$Primary.Cytogenetic.Code, var.axes=FALSE) + ggtitle("Cytogenetic Group")
jpeg("pca_multiplot.jpg", 12, 4, units="in", res=400)
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

# screeplot
jpeg("pca_compimp.jpg", 6, 4, units="in", res=400)
# screeplot(data_pca, main="Variance Explained, by Component", xlab="Component")
fviz_screeplot(res.pca, ncp=10)
dev.off()

# feature/variable contributions
p1 <- fviz_contrib(res.pca, choice = "var", axes = c(1), top=50)
p2 <- fviz_contrib(res.pca, choice = "var", axes = c(1,2), top=50)
p3 <- fviz_contrib(res.pca, choice = "var", axes = c(1,2,3), top=50)
jpeg("pcavar_multiplot.jpg", 12, 12, units="in", res=400)
grid.arrange(p1, p2, p3, ncol=1)
dev.off()

p1 <- fviz_contrib(res.pca, choice = "var", axes = c(1), top=50)
p2 <- fviz_contrib(res.pca, choice = "var", axes = c(2), top=50)
p3 <- fviz_contrib(res.pca, choice = "var", axes = c(3), top=50)
jpeg("pcavar_multiplot_axes123.jpg", 12, 12, units="in", res=400)
grid.arrange(p1, p2, p3, ncol=1)
dev.off()

# pca summary multiplot
p1 <- fviz_screeplot(res.pca, ncp=20)
p2 <- ggbiplot(data_pca, ellipse=TRUE, groups=degfilt.se$Risk.group, var.axes=FALSE) + ggtitle("Biplot (Group = Risk Group)")
p3 <- ggbiplot(data_pca, ellipse=TRUE, groups=degfilt.se$deg.risk, var.axes=FALSE) + ggtitle("Biplot (Group = Binary Risk Group)")
p4 <- ggbiplot(data_pca, ellipse=TRUE, groups=degfilt.se$Primary.Cytogenetic.Code, var.axes=FALSE) + ggtitle("Biplot (Group = Cytogenetic Group)")

jpeg("pcasummary_multiplot.jpg", 12, 10, units="in", res=400)
grid.arrange(p1,p2,p3,p4,layout_matrix=matrix(c(1,1,1,2,3,4), nrow=2, byrow=TRUE))
dev.off()

#---------------
# TSNE Analysis
#---------------
# t-SNE analysis
# Note: from docstring...
# 'perplexity should not be bigger than 3 * perplexity < nrow(X) - 1'
# 'This value effectively controls how many nearest neighbours are taken into account when constructing the embedding in the low-dimensional space.'

# Rtsne(t(data), dims=2, perplexity=50)
tsne <- Rtsne(t(data), perplexity=30)
jpeg("tsne_pxty30_biplot.jpg", 7, 5, units="in", res=400)
plot(tsne$Y, col=ifelse(degfilt.se$deg.risk==1, "red", "blue"),
     main="t-SNE Results")
legend("topright",title="BRG", legend=c("1","0"), col=c("red","blue"), pch=c(16,16), horiz = T)
dev.off()

#-----------------------------------------
# Correlation of Feature Ranked Importance
#-----------------------------------------
# rank correlations
# use actual values for ranks
lasso1rank = rank(st$lasso_coef_rep1)
svm1rank = rank(st$svm1_weights)
xg5rank = rank(st$xg5_imp)
rf10krank = rank(st$rfnb_10k_MeanDecNodeImp)
# append pca ranks
pcdf <- pccont
pcdf <- pcdf[order(match(pcdf$name, st$X)),]
identical(pcdf$name, st$X)
pcvc1 <- rank(pcdf$cont_pca1)
pcvc123 <- rank(pcdf$cont_pca123)
# rankdf
rankdf.all = data.frame(lasso1=lasso1rank,
                     svm1=svm1rank,
                     xg5=xg5rank,
                     rf10k=rf10krank,
                     pcc1=pcvc1,
                     pcc123=pcvc123)

cormat.all <- round(cor(rankdf.all,method="spearman"),3)

# correlation analysis among gene subsets
# among lasso1 genes, with reranking
stsub <- st[!st$lasso_coef_rep1==0,]
pcdfsub <- pcdf[!st$lasso_coef_rep1==0,]

lasso1rank = rank(stsub$lasso_coef_rep1)
svm1rank = rank(stsub$svm1_weights)
xg5rank = rank(stsub$xg5_imp)
rf10krank = rank(stsub$rfnb_10k_MeanDecNodeImp)
pcvc1sub <- rank(pcdfsub$cont_pca1)
pcvc123sub <- rank(pcdfsub$cont_pca123)

rankdf.sub = data.frame(lasso1=lasso1rank,
                        svm1=svm1rank,
                        xg5=xg5rank,
                        rf10k=rf10krank,
                        pcc1=pcvc1sub,
                        pcc123=pcvc123sub)

cormat.lasso1 <- round(cor(rankdf.sub, method="spearman"),3)

# among xg5 genes, with reranking
stsub <- st[!st$xg5_imp==0,]
pcdfsub <- pcdf[!st$xg5_imp==0,]

lasso1rank = rank(stsub$lasso_coef_rep1)
svm1rank = rank(stsub$svm1_weights)
xg5rank = rank(stsub$xg5_imp)
rf10krank = rank(stsub$rfnb_10k_MeanDecNodeImp)
pcvc1sub <- rank(pcdfsub$cont_pca1)
pcvc123sub <- rank(pcdfsub$cont_pca123)

rankdf.sub = data.frame(lasso1=lasso1rank,
                        svm1=svm1rank,
                        xg5=xg5rank,
                        rf10k=rf10krank,
                        pcc1=pcvc1sub,
                        pcc123=pcvc123sub)

cormat.xg5genes <- round(cor(rankdf.sub, method="spearman"),3)

#-----------------------------
# Correlations Visualizations
#-----------------------------
cp1 <- ggcorrplot(cormat.all, hc.order = FALSE, type = "lower",
                  outline.col = "white",
                  legend.title="Rho", title="All DEGs") +
  geom_tile(colour = "black") +
  geom_text(aes(label = round(value, 1)))

cp2 <- ggcorrplot(cormat.lasso1, hc.order = FALSE, type = "lower",
                  outline.col = "white",
                  legend.title="Rho", title="Lasso1 Genes") +
  geom_tile(colour = "black") +
  geom_text(aes(label = round(value, 1)))

cp3 <- ggcorrplot(cormat.xg5genes, hc.order = FALSE, type = "lower",
                  outline.col = "white",
                  legend.title="Rho", title="XGB5 Genes") +
  geom_tile(colour = "black") +
  geom_text(aes(label = round(value, 1)))

# composite plot
jpeg("corhmmulti_consensusimprank.jpg", 12, 4, units="in", res=400)
grid.arrange(cp1,cp2,cp3,layout_matrix=matrix(c(1,2,3), nrow=1, byrow=TRUE))
dev.off()

#---------------------------------
# Upset Plot of "Important" Genes
#---------------------------------
usig <- list()
usig[["lasso1"]] = as.character(st[!st$lasso_coef_rep1==0,]$X)
usig[["svm1"]] = as.character(st[abs(st$svm1_weights)>=0.008,]$X)
usig[["rf10"]] = as.character(st[abs(st$rfnb_10k_MeanDecNodeImp)>=0.1,]$X)
usig[["xg5"]] = as.character(st[!st$xg5_imp==0,]$X)
usig[["pc1"]] = pccont[pccont$above_theoimp_pc1=="Y",]$name
usig[["pc123"]] = pccont[pccont$above_theoimp_pc1=="Y",]$name

for(i in 1:length(usig)){
  names(usig)[i] <- paste0(names(usig)[i]," (",length(usig[[i]])," genes)",collapse="")
}

# upset for filtered gene lists
jpeg("upset_cmlmodelgenes.jpg",7,5,units="in", res=400)
upset(fromList(usig), order.by = "freq", nsets=length(usig),
      sets.bar.color = "blue", main.bar.color = "red")
grid.text("Overlap in Selected and Filtered\nModel Genes", x=0.65, y=0.95, gp=gpar(fontsize=9))
dev.off()

#---------------------
# Compare GSE Results
#---------------------
# get filtered gene id lists
lasso1ids = as.character(st[!st$lasso_coef_rep1==0,]$ensembl_gene_id)
svm1ids = as.character(st[abs(st$svm1_weights)>=0.008,]$ensembl_gene_id)
rf10kids = as.character(st[abs(st$rfnb_10k_MeanDecNodeImp)>=0.1,]$ensembl_gene_id)
xg5ids = as.character(st[!st$xg5_imp==0,]$ensembl_gene_id)

lx = list(lasso1ids, svm1ids, rf10kids, xg5ids)
lplot = c()
lnames = c("lasso1","svm1","rf10k","xg5")
for(i in 1:4){
  lplot[[paste0(lnames[i]," (",length(lx[[i]])," genes)")]] <- lx[[i]]
}

# do GSE comparisons by model type
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
lasso1gse <- enrichr(c(as.character(st[!st$lasso_coef_rep1==0,]$hgnc_symbol)), dbs)
# nrow(st[abs(st$svm1_weights)>=0.008,]) # 81
svm1gse <- enrichr(c(as.character(st[abs(st$svm1_weights)>=0.008,]$hgnc_symbol)), dbs)
# nrow(st[st$rfnb_10k_MeanDecNodeImp>=0.1,]) # 91
rf10kgse <- enrichr(c(as.character(st[abs(st$rfnb_10k_MeanDecNodeImp)>=0.1,]$hgnc_symbol)), dbs)
# nrow(st[!st$xg5_imp==0,]) # 13
xg5gse <- enrichr(c(as.character(st[!st$xg5_imp==0,]$hgnc_symbol)), dbs)

gselist = list("lasso1"=lasso1gse,
               "svm1gse"=svm1gse,
               "rf10kgse"=rf10kgse,
               "xg5gse"=xg5gse)

# upset plots for GSE term overlap
# include filter for gse unadj pval <= 0.05
for(i in 1:3){
  lplot = list()
  whichdf <- i
  fn = c("gomolfun2015","gocellc2015","gobiolproc2015")
  fnplot <- fn[whichdf]
  for(l in 1:length(gselist)){
    namei <- names(gselist)[l]
    dfi <- gselist[[l]][[whichdf]]
    dfi <- dfi[dfi$P.value<=0.05,]
    dati <- dfi$Term
    lplot[[paste0(namei," (",length(dati)," terms)")]] <- dati
    uplotmain = names(gselist[[l]])[whichdf]
  }
  jpeg(paste0("upset_cmlgse_",fnplot,".jpg",collapse=""),7,5,units="in", res=400)
  upset(fromList(lplot), order.by = "freq", nsets=length(lplot),
        sets.bar.color = "blue", main.bar.color = "red")
  grid.text(paste0("GSE Term Overlap in\n", uplotmain), x=0.65, y=0.95, gp=gpar(fontsize=12))
  dev.off()
}

#----------------------------
# GSE Analaysis for All DEGS
#----------------------------
# padj cutoff < 0.01
deggse = enrichr(c(as.character(rowData(degfilt.se)$hgnc_symbol)), dbs)

# filter adj pvals
tnames = c("deggse_gomolfun2015.csv","deggse_gocellc2015.csv","deggse_gobiolproc2015.csv")
for(i in 1:3){
  dfi <- deggse[[i]]
  dfi <- dfi[dfi$Adjusted.P.value<=0.01,]
  deggse[[names(deggse)[i]]] <- dfi
  write.csv(deggse[[i]],file=tnames[i])
}



