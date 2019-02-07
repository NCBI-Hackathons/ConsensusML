ConsensusML [![DOI](https://zenodo.org/badge/169003170.svg)](https://zenodo.org/badge/latestdoi/169003170)
---

# Detecting Cancer Biomarkers from RNA-seq Using Machine Learning
Machine Learning to Detect Cancer Biomarkers from RNAseq Data

- Workflow to apply machine learning methods for feature selection and selection consensus, to determine sets of best discriminating gene biomarkers using RNA-seq data from heterogeneous cancer populations (e.g. pediatric AML, etc.).

- This would be a workflow that is publicly available to bioinformaticians who use RNA-seq gene expression data to characterize tumors, and would assist in finding biomakers that are best able to identify specific subpopulations in a cancer cohort (e.g. TARGET AML, TCGA AML, TARGET NBL, TARGET WT). 

# Background
Leukemia is a cancer of the blood arising in white blood cells of the bone marrow. It poses a substantial population burden as the most common pediatric cancer ([Steliarova-Foucher et al 2017](https://www.ncbi.nlm.nih.gov/pubmed/28410997); [SEER data](https://seer.cancer.gov/statfacts/html/amyl.html)). Acute Myelogenous Leukemeia (AML) is a type of leukemia impacting the myeloblast stem cells. AML arises at a current rate of approximately 20,000 cases per year, with 27.4% 5-year survival 2. It is a molecularly heterogeneous cancer, with several clinically relevant subtypes, including perhaps dozens of subtypes defined by factors ranging from cell differentiation state to cytogenic and sequencing assays ([Yi G. et al 2019](https://www.sciencedirect.com/science/article/pii/S2211124718320588?via%3Dihub); [Tyner J. W. et al 2018](https://www.nature.com/articles/s41586-018-0623-z)). Pediatric AML is characterized at a molecular level by rare somatic mutations, absence of common adult AML mutations, and relatively frequent structural variants ([Bolouri H et al 2018](https://www.nature.com/articles/nm.4439)). Here, we apply several machine learning approaches for feature selection of RNA-seq data from both pediatric and adult AML cases. Our goal was to better understand gene expression-based heterogeneity underlying AML cases, as well as age-related and -unrelated dysregulation patterns. We used clinical and assay data from pediatric cancer patients from the Therapeutically Applicable Research To Generate Effective Treatments (TARGET) initiative (https://ocg.cancer.gov/programs/target/).

# Methods and Analysis Overview

We were interested in applying machine learning principles for feature selection, to identify the most important genes and gene sets for predicting clinically-relevant classifiers in pediatric and adult AML cases. Classifiers of main interest include age, stage, and survival. For this investigation, we focused on the risk group sample classifier. We performed both pan-cancer and cancer-specific analyses of TARGET pediatric cancers. For analysis of AML cases, we combined primary peripheral blood and bone marrow samples. 

## Methods Workflow
!["ML RNA-seq biomarkers, methods flow chart"](https://github.com/NCBI-Hackathons/RNAseq_Cancer_Biomarkers/blob/master/methods.jpg "Day 1 Flowchart")

## Expression Data
TARGET gene counts were obtained from the Genomic Data Commons (https://gdc.cancer.gov/), which are based on RNA-seq run using the Illumina HiSeq platform. Gene counts were normalized using trimmed mean of m-values (TMM) method. We further pre-filtered TMM-normalized expression based on extent of differential expression between these classifiers of interest, using multiple thresholds. After these preprocessing and QC methods were complete, for each classifier of interest we randomly divided data in training and test subsets, conserving classifier frequency in each subset.

## Analysis Approach using Machine Learning Algorithms
We then applied an "ensemble" learning approach comparing feature selection results. We assessed results using multiple machine learning algorithms from various R and Python packages, including: 1. Support Vector Machines (SVM) using [e1071](https://cran.r-project.org/web/packages/e1071/index.html); 2. Random Forest with boosting in Python; 3. Neural Networks with [keras](https://cran.r-project.org/web/packages/keras/index.html); 4. Logistic regression in R; 5. Elastic net with Lasso using [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html); 6. [AutoML](https://pypi.org/project/automl/). 

# Results and Feature Selection

Using normalized and pre-filtered RNA-seq expression data, we fitted models using various algorithms as described in Methods. For models that performed well on the filtered gene set, we then identified the most important gene features for model prediction. We assessed consensus of selected features across algorithm classes, and for the most common recurrently selected features, we mined the scientific literature for evidence validating these genes' functional roles in leukemia and AML.

# Links to Shared Documents

## 1. Manuscript
https://docs.google.com/document/d/1DPAmUFfggAnAjsMIPTs1hV90k25ZKckYLi18b3dBot0/edit#

## 2. Day 2 Presentation
https://docs.google.com/presentation/d/1HxHyaGLNxAbhsEd2OVs6R3HiGGu5hrJO_lcLS5ffZlc/edit#slide=id.g4f487fb995_0_278

## 3. Final Hackathon Presentation
https://docs.google.com/presentation/d/1WdvgktxKQAXjABnFVWuMOW6ENCd5gwZKb_SjKrnLZbs/edit#slide=id.g4ec347d104_0_7
