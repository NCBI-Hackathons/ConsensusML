**Project:** ConsensusML

**Author:** Jenny Smith


The code in this directory provides the methods for
1. Downloading Gene Expression Data from the Genomic Data Commons
  - GDC_Data_Download.Rmd
2. Normalization and Batch Effect Investigations
  - Normalization_and_Batch_Effect_Investigation.Rmd
3. Differential expression and LASSO logistic regression using differentially expressed genes (DEGs) as variables
  - Differential_Expression_and_Lasso.Rmd

The focus of these Investigations are to identify genes that are associated with high and standard risk AML, first by differential expression analysis and then to further select genes from a large list of DEGs using the LASSO regression.  These would be genes that are associated with high and standard risk clinical features.  


This procedure produced a small  number of genes which can be easily further investigated for defining theraputic targets/biomarkers or potentially prediction of poor prognosis using diagnostic samples.

!["volcano_plot"](https://github.com/NCBI-Hackathons/RNAseq_Cancer_Biomarkers/blob/master/JSmith_code/Results/RiskGroup_DEGs_VolcancoPlot.tiff "volcano_plot")

!["ML RNA-seq biomarkers, methods flow chart"](https://github.com/NCBI-Hackathons/RNAseq_Cancer_Biomarkers/blob/master/methods.jpg "Day 1 Flowchart")
