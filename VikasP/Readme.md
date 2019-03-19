**Implementation using MLseq and DEseq R packages**

script: AML.R

files required to run: TARGET_NBL_AML_RT_WT_HTSeq_Counts.csv.gz, TARGET_AML_High.Std.Risk_vs_LowRisk_DEGs.csv, TARGET_AML_High.Std.Risk_vs_LowRisk_DEGs.csv

Detailed instructions for using MLseq and DEseq as used here: https://www.bioconductor.org/packages/release/bioc/vignettes/MLSeq/inst/doc/MLSeq.pdf

My approach was to use the MLseq and DEseq R packages from bioconductor. Data from TARGET_NBL_AML_RT_WT_HTSeq_Counts.csv.gz contains expression data for the full ~60k genes. This file was subsetted to only include samples with clinical data from TARGET_AML_High.Std.Risk_vs_LowRisk_DEGs.csv. The models used were voomNSC, voomDLDA, PLDA, PLDA2, and nblda. Results from this are available in VikasP as plda.txt, plda2.txt, and voomNSC.txt. 

Afterwards TARGET_NBL_AML_RT_WT_HTSeq_Counts.csv.gz was subsetted again using only genes available in, TARGET_AML_High.Std.Risk_vs_LowRisk_DEGs.csv, which were selected as demonstrating differential gene expression. The same models as above were retrained on this. Data is available in VikasP as DEG_plda.txt, DEG_plda2.txt, and DEG_voomNSC.txt. 

Most of the genes identified by the PLDA 2 model are related to various anemias, alpha/beta Thalassemias, Hemoglobin C/D/E/M. Interestingly CSF3R, associated with chronic myeloid leukemia, and GATA2 associated with Acute Myeloid Lukemia, were found with the PDLA/PDLA2 models. 
