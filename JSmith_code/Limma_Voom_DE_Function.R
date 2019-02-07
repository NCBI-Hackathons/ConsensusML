#Jenny Smith

#April 5, 2017

#purpose: to perform differential expression using limma voom transformed data



voom_DE <- function(counts.df, ref, pheno){
  #counts.df is a dataframe with count data, with genes as rownames
  #pheno is a character vector with patient IDs as names, and the status for each in each group(eg pos,neg)
  require(edgeR)
  library(limma)
  
  #ensure correct order for both expn and counts.df
  samples <- intersect(names(pheno), colnames(counts.df))
  pheno <- pheno[samples]
  counts.df <- counts.df[,samples]
  
  
  groups <- unique(pheno)
  groups <- c(groups[groups != ref], ref) #order so that reference is second 
  pheno.f <- factor(pheno, levels=groups)

  dge <- DGEList(counts = counts.df, group = pheno.f)

  keep.dge <- rowSums(cpm(dge) >= 1) > (0.05*ncol(counts.df)) #5% of samples with CPM >= 1
  dge <- dge[keep.dge,]
  dge <- calcNormFactors(dge)

  design <- model.matrix(~0 + pheno.f, data=dge$samples)
  colnames(design) <- levels(pheno.f)
  cont.matrix <- makeContrasts(contrasts = paste(groups, collapse = "-"), levels = design)
  
  
  v.lv <- voom(dge, design, plot = FALSE)
  

  fit <- lmFit(v.lv, design)
  fit <- contrasts.fit(fit, contrasts = cont.matrix)
  fit <- eBayes(fit)
  table <- topTable(fit, number = 20000, p.value=0.05, adjust.method="BH", sort.by="P",lfc=1)
  


  list <- list(design, v.lv, fit, table)
  names(list) <- c("desingMatrix", "voomTransformation", "fit", "DEGs")
  return(list)
}



#Limma Voom
# When the library sizes are quite variable between samples, then the voom approach is theoretically
# more powerful than limma-trend.https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pcounts.df















